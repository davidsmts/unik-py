from Bio import SeqIO
import numpy as np


def get_reads_from(seq, read_len=200, cov=5):
    L = len(seq)
    read_amt = int(L/read_len * cov)
    reads = []
    startposis = []
    for _ in range(read_amt):
        start = np.random.randint(0,L-read_len)
        read = seq[start:start+read_len]
        reads.append(read)
        startposis.append(start)
    return reads, startposis


def get_skmer(sequence, profile):
    spaced_kmer = sequence * profile
    spaced_kmer = spaced_kmer[spaced_kmer != 0]
    s = int(''.join(str(x) for x in spaced_kmer))
    return s


def parse_nucleotides(sequence):
    new_seq = []
    map_to_vals_L = {"A": 1, "C": 2, "G": 3, "T":4}
    map_to_vals_S = {"a": 1, "c": 2, "g": 3, "t":4}
    for symbol in sequence:
        if symbol.isupper():
            new_seq.append(map_to_vals_L[symbol])
        else:
            new_seq.append(map_to_vals_S[symbol])
    
    return new_seq


def read_fasta(filename):
    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    reads = []
    headers = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        reads.append(sequence)
        headers.append(name)
    return reads, headers


def save_clusters(clusters, all_reads):
    clustered_reads = {}
    for (i, el) in zip(range(len(clusters)), clusters):
        if el not in clustered_reads:
            clustered_reads[el] = [el]
        else:
            clustered_reads[el].append(i)
    files = []
    for key in clustered_reads.keys():
        read_ids = clustered_reads[key]
        filename = "../data/clusters_CAMI2/"+str(key)+".fasta"
        files.append(filename)
        with open(filename, "w") as file:
            for ID in read_ids:
                file.write(">"+str(ID)+"\n")
                file.write("".join(all_reads[ID]) + "\n")      