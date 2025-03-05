import argparse
import unik

parser = argparse.ArgumentParser(description='Python version of Unik metagenome assembler')

# Required positional arguments
parser.add_argument('files', type=list, help='List of files.')
parser.add_argument('target_dir', type=int, help='The directory where the clusters will be written.')

# Optional arguments
#parser.add_argument('--kmer-profile', type=str, help='The spaced k-mer profile. Is 1111110110110101110101011101011111101111 by default.')

# Boolean switch
#parser.add_argument('--has-reads', action='store_true', help='Flag indicating wether the fasta files contain reads or whole genomes.')

args = parser.parse_args()
unik.unik(args.files, args.target_dir)