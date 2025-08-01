from optparse import OptionParser
from Readfile import *
from Edmonds import *
import os,sys
import subprocess
#get the absolute path of input file
def getPath(path):
    path1=path.split("/")
    if path1[0] == ".":
        if (len(path1)==1):
            newpath=os.getcwd()
        else:
            newpath=os.getcwd()
            for ele in path1:
                if ele !=".":
                    newpath=newpath+"/"+ele
    elif path1[0]=="..":
        i = 0
        for ele in path1:
            if ele == "..":
                i=i+1
        path2=os.getcwd()
        path2=path2.split("/")
        newpath="/"+path2[0]
        if len(path2)-i > 1:
            for j in range(1,len(path2)-i):
                newpath=newpath+"/"+path2[j]
        for j in range(i,len(path1)):
            newpath=newpath+"/"+path1[j]
    else:
        newpath=path
    return newpath


def main():
    usage = "usage: python %prog <-P path> <-I input> <-D datatype>"
    description = "Input integer copy number profile. Columns correspond to chromosomal position. Rows correspond to cells."
    op = OptionParser(version="%prog 1.0",description=description,usage=usage,add_help_option=False)
    op.add_option("-h","--help",action="help",
                  help="Show this help message and exit.")
    op.add_option("-P","--Path",dest="Path",type="str",
                  help="Path to script")
    op.add_option("-I","--Input",dest="Input",type="str",
                  help="Input file")
    op.add_option("-G","--Genome",dest="Genome",type="str",
                  help="Genome version hg19 or hg38")
    op.add_option("-O","--Output",dest="Output",type="str",
                  help="Output path")
    op.add_option("-D","--Datatype",dest="Datatype",type="str",
                  help="The type of input data. Either D (DNA-seq) or R (RNA-seq).")
    op.add_option("-W","--Windows",dest="Windows",type="str",
                  help="the number of genes you want to merge when you input copy number profile inferred from scRNA-seq. Default 30.")
    op.add_option("-R","--Permutation",dest="Permutation",type="str",
                  help="Whether reconstructed permuted tree (T) or not (F). If not, permuted copy number profile will be used to perform LSA. Default value is F due to time cost.")

    (options,args) = op.parse_args()
    # check input parameters. Package path, input file, data type and genome version are required.
    if not options.Path or not options.Input or not options.Datatype or not options.Genome:
        op.print_help()
        sys.exit(1)

    # get the input parameters
    currentPath=os.getcwd()
    scTreepath=options.Path
    scTreepath=getPath(scTreepath)
    filename=options.Input
    filename=getPath(filename)
    hg=options.Genome
    if not options.Output:
        outpath=currentPath
    else:
        outpath=options.Output
        outpath=getPath(outpath)
    os.system("mkdir -p "+outpath)
    datatype=options.Datatype
    os.system("mkdir -p "+outpath+"/temp")
    writename=outpath+"/CNV.tree.txt"
    
    if not options.Permutation:
        permutation = "F"
    else:
        permutation=options.Permutation
    print "Transfer data to segmental level"

    #change the input copy number profile to the matrix format used to infer Tree
    #if the data type = R, estimate integer copy number profile by merging the number of adjacent genes.
    #Default number of genes is 30.
    # Run R script from the original directory to avoid path issues
    if datatype == "D":
        os.system("Rscript "+scTreepath+"dataTransfer.R "+filename+" "+datatype)
    elif datatype == "R":
        if not options.Windows:
            print "The number of genes which are merger into the bin is default value 30. If you want change it please specify the value through -W"
            delt = str(30)
        else:
            delt=options.Windows
        os.system("Rscript "+scTreepath+"dataTransfer.R "+filename+" "+datatype+" "+scTreepath+" "+delt)
    else:
        print "Please provide the correct inputfile type through -D either 'D' or 'R'."
    CNVfile=filename+".CNV.txt"
    
    # Now change to temp directory for the rest of the processing
    os.chdir(outpath+"/temp")
    # Copy the CNV file to the temp directory so it can be found
    # Use absolute path to the processed CNV file
    abs_CNVfile = os.path.join(currentPath, os.path.basename(CNVfile))
    print(os.getcwd())
    os.system("cp ../../"+filename+".CNV.txt"+" .")
    print(os.listdir('.'))
    print "Inferring MEDALT."

    #Identifying root node from input file.
    #If a diploidy genome is not input, will add an extra diploidy node as root
    # Use just the filename without path since we're in the temp directory
    CNVfile_basename = os.path.basename(CNVfile)
    (nodes,root) = read(CNVfile_basename)
    node_name_list = nodes.keys()
    print 'Read file.'
    #calculation of MED distance
    from ComputeDistance import dist
    g = create_tree(nodes, node_name_list,root)
    print 'Tree crated.'
    #Inference of tree and outp
    from mdmst_original import compute_rdmst as compute_rdmst_original
    # Replace with new tree
    from mdmst import compute_rdmst
    result = compute_rdmst(g, root)
    print 'RDMST computed'
    # Use basename for the write file since we're in the temp directory
    write_basename = os.path.basename(writename)
    write=open(write_basename,'w')
    tree=result[0]
    out1="from"+"\t"+"to"+"\t"+"dist"
    print >> write, out1
    for ele in tree.keys():
        out=ele
        for value in tree[ele].keys():
            out1=out+"\t"+value+"\t"+str(tree[ele][value])
            print >> write,out1
    write.close()
    # Copy the output file back to the output directory
    print(write_basename)
    print(writename)
    os.system("cp CNV.tree.txt ..")
    print "MEDALT inferrence finish."

    #Permutation process for lineage speciation analysis (LSA)
    #Permutate the copy number profile by chromosome into different cells
    #if permutation == T, tree corresonds to each permuted data will be inferred by above algorithm
        # we infer 100 permutation tree due to the cost of time
    #if permutation == F, just the permutation datasets rather than permutation tree are used to estimate the significance.
        # we permute 500 times.
    #This process is done by R code.
    if permutation == "T":
        permutationPath=outpath+"/temp"
        print "Reconstructing tree based on permutation data."
        print "This will take a long time! Please have some coffee."

        #permute copy number profile
        if datatype == "D":
            os.system("Rscript "+scTreepath+"permutationCNA.R "+scTreepath+" "+filename+" "+datatype+" "+permutationPath)
        elif datatype == "R":
            os.system("Rscript "+scTreepath+"permutationCNA.R "+scTreepath+" "+filename+" "+datatype+" "+permutationPath+" "+delt+" "+hg)

        #Infer permutation tree
        for j in range(1,101):
            permutefile=permutationPath+"/permute."+str(j)+".CNV.txt"
            (nodes,root) = read(permutefile)
            node_name_list = nodes.keys()
            g = create_tree(nodes, node_name_list,root)
            result = compute_rdmst(g, root)
            permuteTree=permutefile+".celltree.txt"
            write=open(permuteTree,'w')
            tree=result[0]
            out1="from"+"\t"+"to"+"\t"+"dist"
            print >> write, out1
            for ele in tree.keys():
                out=ele
                for value in tree[ele].keys():
                    out1=out+"\t"+value+"\t"+str(tree[ele][value])
                    print >> write,out1
            write.close()
        print "Pemutation tree finish."
        print "Performing LSA."

        #Identifying CNAs associated with cellular lineage expansion.
        command = "Rscript "+scTreepath+"LSA.tree.R "+scTreepath+" "+filename+" "+writename+" "+CNVfile+" "+outpath+" "+datatype+" "+hg+" "+permutationPath
        os.system(command)
    elif permutation == "F":
        print "Performing LSA."
        # Append filename and CNVfile to the input folder from args. Add sctreepath and the foldername of the input filename
        filename_folder = os.path.dirname(filename)
        output_folder = os.path.dirname(writename)
        CNVfile = os.path.join(scTreepath, filename_folder, os.path.basename(CNVfile))
        filename = os.path.join(scTreepath, filename_folder, os.path.basename(filename))
        writename = os.path.join(scTreepath, output_folder, os.path.basename(writename))
        print('CNVfile is:', CNVfile)
        print('filename is:', filename)
        print('writename is:', writename)
        command="Rscript "+scTreepath+"LSA.tree.R "+scTreepath+" "+filename+" "+writename+" "+CNVfile+" "+outpath+" "+datatype+" "+hg
        print(command)
        os.system(command)
    print "Done!"

if __name__ == "__main__":

    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me, see you!\n")
        sys.exit(0)
        sys.exit(0)
