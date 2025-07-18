#main R function to perform linegae speciation analysis
#input: 1. the path of R script LSAfunction.R
#       2. input copy number profile from scDNA-seq or scRNA-seq
#       3. the infer tree file based on integer copy number profile
#       4. integer copy number matrix
#       5. output path of results
#       6. datatype "D" (DNA-seq) or "R" (RNA-seq)
#       7. genome version: hg19 or hg38
#       8. optional. If you have the permutation tree, input the path


#Check the installation of dependency packages: HelloRangers and igraph
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
if (!"HelloRanges" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("HelloRanges")
}
if (!"igraph" %in% installed.packages()){
  install.packages("igraph")
}
library(HelloRanges)
library(igraph)

#read input
args<-commandArgs(T)
datapath=args[1]
inputfile=args[2]
treeName=args[3]
CNVfile=args[4]
outpath = args[5]
datatype=args[6]
hg=args[7]
if (length(args)==8){
  permutationPath=args[8]
}
source(paste(datapath,"/LSAfunction.R",sep=""))
set.seed(1234)

#Plot cell tree figure
print.noquote("Visualization MEDALT!")
celltree=read.csv(treeName,sep="\t")
nodes=data.frame(id=union(as.character(celltree[,1]),as.character(celltree[,2])),size=5)
nodes$color="lightblue"
nodes$color[nodes$id==setdiff(as.character(celltree[,1]),as.character(celltree[,2]))]="black"
net <- graph_from_data_frame(d=celltree, vertices=nodes, directed=T)
pdf(file="singlecell.tree.pdf",width = 5,height = 5,useDingbats = F)
plot(net, vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label=NA)
dev.off()

#read input copy number data
data = read.csv(inputfile,sep="\t",header = TRUE)

#read reference genome file
if (hg=="hg19"){
  reference=read.csv(paste(datapath,"/gencode_v19_gene_pos.txt",sep=""),sep="\t",header=F)
  refer.band=read.csv(paste(datapath,"/hg19.band.bed",sep=""),sep="\t",header = F)
}
if (hg=="hg38"){
  reference=read.csv(paste(datapath,"/gencode_v38_gene_pos.txt",sep=""),sep="\t",header=F)
  refer.band=read.csv(paste(datapath,"/hg38.band.bed",sep=""),sep="\t",header = F)
}

#generate the genome position data at arm and band level
chrom=as.character(unique(refer.band[,1]))
arm.band=c()
band.region=c()
for (i in chrom){
  subrefer=refer.band[refer.band[,1]==i,]
  parm=subrefer[grep("p",subrefer[,4]),]#p arm region
  qarm=subrefer[grep("q",subrefer[,4]),]#q arm region
  armregion=data.frame(chr=i,start=c(parm[1,2],qarm[1,2]),end=c(parm[dim(parm)[1],3],qarm[dim(qarm)[1],3]),band=c("p","q"))
  arm.band=rbind(arm.band,armregion)
  subband=do.call(rbind,strsplit(as.character(subrefer[,4]),split="[.]"))
  bandname=unique(subband[,1])
  pos=sapply(bandname,function(x,subband,subrefer){
    k=which(subband[,1]==x)
    start=subrefer[min(k),2]
    end=subrefer[max(k),3]
    return(c(start,end))
  },subband,subrefer)
  subregion=data.frame(chr=i,start=pos[1,],end=pos[2,],band=bandname)
  band.region=rbind(band.region,subregion)
}
arm.band$length=arm.band$end-arm.band$start
band.region$length=band.region$end-band.region$start
refer.band$ID=paste(as.character(refer.band$V1),as.character(refer.band$V4),sep=":")

#generate the bin position based on input copy number data
print.noquote("LSA segmentation!")
if (datatype=="D"){
  region=data[,1:2]
  region[,1]=paste("chr",region[,1],sep="")
  region$end=region[,2]
  colnames(region)=c("chrom","chrompos","end")
}else if (datatype=="R"){
  data=round(data*2)
  index=match(row.names(data),as.character(reference[,1]))
  newdata=cbind(reference[index[!is.na(index)],2:3],data[!is.na(index),])
  rownames(newdata)=rownames(data)[!is.na(index)]
  newdata=newdata[as.character(newdata[,1])!="chrM"&as.character(newdata[,1])!="chrY",]
  newdata[,1]=as.character(newdata[,1])
  newdata[newdata[,1]=="chrX",1]="chr23"
  region=newdata[,1:2]
  region[,3]=region[,2]
  colnames(region)=c("chrom","chrompos","end")
  cnv=t(data)
  data=newdata
}
write.table(region,"region.bed",col.names = F,row.names = F,sep="\t",quote = F)

##correspond input genomic bin to genome band ID
code=bedtools_intersect(paste("-a ",datapath,"/",hg,".band.bed -b region.bed",sep=""))
ans <- eval(code)
ans=as.data.frame(ans)
ans$ID=paste(ans$seqnames,":",ans$name,sep="")
ans$ID1=paste(ans$seqnames,"_",ans$end,sep = "")
region$ID1=paste(region$chrom,"_",region$chrompos,sep="")
ID=unique(ans$ID)

#generate copy number matrix based on genome bandID
newCNV=do.call(cbind,lapply(ID, function(id,data,ans,region){
  chrID=ans$ID1[ans$ID==id]
  index=match(chrID,region$ID1)
  if (length(index)==1){
    return(as.numeric(data[index,3:dim(data)[2]]))
  }else{
    return(round(apply(data[index,3:dim(data)[2]],2,mean)))
  }
},data=data,ans=ans,region))
colnames(newCNV)=ID

#perform lineage speciation analysis
print.noquote("Calculating CFL")

#calculate the depth, the number of children nodes for each cell
cell=union(as.character(celltree[,1]),as.character(celltree[,2]))
cell=data.frame(cell=cell)
cell$depth=sapply(as.character(cell$cell),depthFunction,cellTree=celltree)
cell$subtreesize=sapply(as.character(cell$cell),subtreeSize,cellTree=celltree)

#perform lineage speciation analysis for cells that the number of children is no less than 5.
cell1=cell[cell$subtreesize>=5,]
cell1=cell1[cell1$cell!="root",]

#calculte CFL at genomic bin level
Gscore=lapply(as.character(cell1$cell),lineageScore,newCNV,celltree)
names(Gscore)=as.character(cell1$cell)

#read genes from 11 oncogenic pathway for estimation of gene level
pathwaygene=read.csv(paste(datapath,"/pathwaygene.txt",sep=""),sep="\t")
index=match(pathwaygene$name,reference[,1])
pathwaygene=data.frame(chr=reference[index[!is.na(index)],2],start=reference[index[!is.na(index)],3],end=reference[index[!is.na(index)],4],name=pathwaygene$name[!is.na(index)],pathway=pathwaygene$pathway[!is.na(index)])

#Gene level copy number profile
if (datatype=="D"){
  cnv=read.csv(CNVfile,sep="\t")
  oncogenicCNV=do.call(cbind,lapply(1:dim(pathwaygene)[1],geneCNAfunction,pathwaygene=pathwaygene,ancestorCNV=cnv,generegion=region))
  colnames(oncogenicCNV)=as.character(pathwaygene$name)
  rownames(oncogenicCNV)=rownames(cnv)
}else if (datatype == "R"){
  data=data[,3:dim(data)[2]]
  index=match(as.character(pathwaygene$name),rownames(data))
  oncogenicCNV=t(data[index[!is.na(index)],])
  oncogenicCNV=round(oncogenicCNV)
  colnames(oncogenicCNV)=as.character(pathwaygene$name)[!is.na(index)]
}
index=apply(oncogenicCNV,2,function(x){
  if (NA %in% x){
    return(0)
  }else{
    return(1)
  }
})
oncogenicCNV=oncogenicCNV[,index==1]

#Calculate CFL at gene level
geneGscore=lapply(as.character(cell1$cell),lineageScore,oncogenicCNV,celltree)
names(geneGscore)=as.character(cell1$cell)
realres=list(cell=cell1,bandGscore=Gscore,geneGscore=geneGscore)

#do calculation for permutation dataset
print.noquote("Calculating permutation CFL")
if (length(args) < 8){
  #if there was no permutation tree, calculate CFL in permutation dataset based on real tree structure
  times=500
  permuteres=lapply(1:times,function(j,data,ID,ans,datatype,pathwaygene,generegion,reference,celltree){
    score=permuteScore(data,ID,ans,datatype,pathwaygene,generegion=region,reference,celltree)
    return(score)
  },data,ID,ans,datatype,pathwaygene,generegion=region,reference,celltree)
}else if (length(args)==8){
  #if there are permutation trees, calculate CFL in permutation dataset based on permuted tree structure
  permutefile=list.files(permutationPath)
  permutefile=permutefile[grep("celltree",permutefile)]
  times=length(permutefile)
  print.noquote(paste("There are ",length(permutefile)," permutation trees."))
  if (length(permutefile)>0){
    permuteres=lapply(permutefile,permuteTreeScore,ID,ans,datatype,pathwaygene,generegion=region,reference,permutationPath)
  }
}

#Estimate emperical p value
print.noquote("Estimate emperical p value")
realcell=realres$cell
realband=realres$bandGscore
realgene=realres$geneGscore
pvalue=lapply(1:dim(realcell)[1],significanceLevel,realband,realgene,permuteres,realcell)
if (length(args) < 8){
  #estimate emperical p value at genomic bin and gene level
  #if there is no permutation tree, default cutoff of pvalue is 0.01.
  bandsig=CollectAsso(pvalue,cutoff=0.01,celltree,realcell)$bandres
  genesig=CollectAsso(pvalue,cutoff=0.01,celltree,realcell)$generes
}else if (length(args)==8){
  #if there are permutation trees, default cutoff of pvalue is 0.05.
  bandsig=CollectAsso(pvalue,cutoff=0.05,celltree,realcell)$bandres
  genesig=CollectAsso(pvalue,cutoff=0.05,celltree,realcell)$generes
}
bandsig=do.call(rbind,lapply(unique(as.character(bandsig$cell)),mergeCNA,bandsig,band.region,arm.band,refer.band))
bandsig=unique(bandsig)
genesig=unique(genesig)
LSAres=list()
if (!is.null(bandsig)){
  index=match(as.character(bandsig$cell),as.character(realcell$cell))
  bandsig$subtreesize=realcell$subtreesize[index]
  cellsig=do.call(rbind,lapply(as.character(unique(bandsig$cell)),CombineRegion,bandsig,refer.band))
  LSAres$bandLSA=cellsig
  paraBand=table(as.character(bandsig$region))
  paraBand=paraBand[paraBand>1]#CNAs of genomic bin associated with more than one independent lineages
  if (length(paraBand)>0){
    paraBandsig=GenePara(bandsig,permuteres,type="band",realcell)#parallel evolution test
    if (!is.null(paraBandsig)){
      LSAres$paraBand=paraBandsig
    }
  }
}
if (!is.null(genesig)){
  index=match(as.character(genesig$cell),as.character(realcell$cell))
  genesig$subtreesize=realcell$subtreesize[index]
  LSAres$geneLSA=genesig
  paraGene=table(as.character(genesig$region))
  paraGene=paraGene[paraGene>1]#gene associated with more than one independent lieage
  if (length(paraGene)>0){
      paraGenesig=GenePara(genesig,permuteres,type="gene",realcell)#parallele evolution estimation
      LSAres$paraGene=paraGenesig
  }
}
print.noquote("Estimate parallel evolution")

#collect all significant results
allsig=c()
if ("geneLSA" %in% names(LSAres)){
  geneLSA=LSAres$geneLSA
  geneLSA$CNA="AMP"
  geneLSA$CNA[geneLSA$Score<0]="DEL"
  allsig=rbind(allsig,geneLSA)
  write.table(geneLSA,paste("gene.LSA.txt",sep=""),col.names = T,row.names=F,sep="\t",quote=F)
}else{
  print.noquote("No LSA is identified at gene level!")
}
if ("bandLSA" %in% names(LSAres)){
  bandLSA=LSAres$bandLSA
  bandLSA$CNA="AMP"
  bandLSA$CNA[bandLSA$Score<0]="DEL"
  allsig=rbind(allsig,bandLSA)
  write.table(bandLSA,paste("segmental.LSA.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
}else{
  print.noquote("No segmental LSA is identified!")
}

#collect all parallel evolution events
paraEvent=c()
if ("paraBand" %in% names(LSAres)){
  paraEvent=rbind(paraEvent,LSAres$paraBand)
}
if ("paraGene" %in% names(LSAres)){
  paraEvent=rbind(paraEvent,LSAres$paraGene)
}
if (!is.null(paraEvent)){
  paraEvent=paraEvent[!is.na(paraEvent$pvalue),]
  if (dim(paraEvent)[1]!=0){
    write.table(paraEvent,paste("parallel.LSA.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=FALSE)
  }
}

#plot LSA Tree
if (!is.null(allsig)){
  LSAnetwork=CNAconnect(allsig,celltree)
  nodes=data.frame(id=union(LSAnetwork[,1],LSAnetwork[,2]),size=5)
  tab=table(as.character(allsig$cell))
  index=match(nodes$id,names(tab))
  nodes$size[!is.na(index)]=nodes$size[!is.na(index)]*tab[index[!is.na(index)]]/5#define node size
  nodes$size[nodes$size<=5]=5
  nodes$color="gray"
  nodes$color[!is.na(index)]=rainbow(length(unique(allsig$cell)))
  annotation=c()
  for (i in 1:dim(nodes)[1]){
    if (as.character(nodes$id[i]) %in% as.character(allsig$cell)){
      CNA=allsig[as.character(allsig$cell)==as.character(nodes$id[i]),]
      CNA=CNA[order(CNA$pvalue),]
      CNA=paste(as.character(CNA$region),as.character(CNA$CNA),sep=":")
      CNA1=CNA[1]
      if (length(CNA)>1){
        for (j in 2:min(3,length(CNA))){
          CNA1=paste(CNA1,CNA[j],sep=";")
          }
      }
      annotation[i]=CNA1
    }
  }
  nodes$annotation=annotation
  nodes$size=nodes$size/max(nodes$size)*30
  links=data.frame(from=LSAnetwork[,1],to=LSAnetwork[,2],weight=as.numeric(LSAnetwork[,3]))
  pdf("LSA.tree.pdf",width = 6,height = 6,useDingbats = F)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
  plot(net, layout=layout_as_tree,vertex.frame.color=NA,vertex.color=nodes$color,edge.arrow.size=.2,vertex.label.cex=0.5,vertex.label=nodes$annotation)
  dev.off()
}
