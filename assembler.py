import subprocess
import sys
from os import listdir
import os.path as op

def spades_assembly_fromclustfiles(directory="../data/clusters_CAMI2/"):
    files = listdir(directory)
    print(files)
    file = files[0]
    #res = subprocess.run(["../plass/build/bin/penguin", "nuclassemble", "--remove-tmp-files", "1", "--num-iterations", "12", "--keep-target", "0", "--contig-output-mode", "1", "../data/clusters/"+file, "../data/assemblies/ass_"+file, "tmp"], capture_output=True)
    print(files)
    for (i, file) in zip(range(len(files)), files):
        if file == ".DS_Store":
            continue
        res = subprocess.run(["../SPAdes-4.1.0-Darwin/bin/spades.py", "-s" , directory+file, "--isolate", "--only-assembler", "-o", directory + "results/" + str(i)+"/"])
    
    print("DONE")


def sga_assembly_fromclustfiles(directory):
    files = listdir(directory)
    print(files)
    for (i, file) in zip(range(len(files)), files):
        if file == ".DS_Store":
            continue
        res0 = subprocess.run(["sga", "assemble" , "../data/clusters_CAMI2/"+file, "-o", "../data/clusters_CAMI2/results_sga/"+str(i)+"/"])
        #...
        res5 = subprocess.run(["sga", "assemble" , "../data/clusters_CAMI2/"+file, "-o", "../data/clusters_CAMI2/results_sga/"+str(i)+"/"])
        


def merge_assembly_results(directory, newfile):
    files = listdir(directory)
    with open(newfile, "w") as nf:
        for file in files:
            if file == ".DS_Store":
                continue
            contig_file = directory+file+"/contigs.fasta"
            if op.isfile(contig_file):
                with open(contig_file, "r") as f:
                    for line in f:
                        nf.write(line)
            else:
                print("nofile " + contig_file)


def eval_assembly_quast(directory):
    newfile = "../data/assemblies/assembly_CAMI2/merged_assembly.fasta"
    merge_assembly_results(directory, newfile)
    res = subprocess.run(["python3", "../../other/quast-5.3.0/quast.py", "../data/assemblies/assembly_CAMI2/merged_assembly.fasta", "-r", "../data/genomes/CAMI2_source_genomes/MIKI-NS02.fasta", "-o", "../data/out/CAMI2/"], capture_output=True) 
    print(res)
