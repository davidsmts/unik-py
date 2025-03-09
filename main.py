import argparse
import unik
import assembler

parser = argparse.ArgumentParser(description='Python version of Unik metagenome assembler')

# Required positional arguments
parser.add_argument('mode', type=str, help="Cluster only, assemble, align or cluster+assemble+align")
parser.add_argument('source', type=str, help='Source directory or file.')
parser.add_argument('target_dir', type=str, help='The directory where the clusters will be written.')

# Optional arguments
#parser.add_argument('--kmer-profile', type=str, help='The spaced k-mer profile. Is 1111110110110101110101011101011111101111 by default.')

# Boolean switch
parser.add_argument('--hasreads', action='store_true', help='Flag indicating wether the fasta files contain reads or whole genomes.')

args = parser.parse_args()

if args.mode == "cluster":
    unik.unik([args.source], args.target_dir, has_reads=args.hasreads)
elif args.mode == "assemble":
    assembler.spades_assembly_fromclustfiles(args.source)
elif args.mode == "align":
    assembler.eval_assembly_quast()
elif args.mode == "all":
    tmp_name = "./tmp/"
    unik.unik(args.files, tmp_name)
    assembler.spades_assembly_fromclustfiles(directory="")
    assembler.eval_assembly_quast()
else: 
    print("Unknown Mode")