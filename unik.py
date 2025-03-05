import numpy as np
import matplotlib.pyplot as plt
import mmh3
from collections import defaultdict, deque
import os.path as op
import misc
import unique_kmers
import count_correction as cc
import clust


def build_hashmaps(all_seqs, all_startposis, profile):
    # BUILD THE HASHMAP
    starts = {}
    ends = {}
    first_appearance = {}
    seqs_kmers = {}
    print("Building hashmap")
    for (index, sequence) in zip(range(len(all_seqs)), all_seqs):
        for i in range(len(sequence) - f +1):
            spaced_kmer = sequence[i:i+f] * profile
            spaced_kmer = spaced_kmer[spaced_kmer != 0]
            s = ''.join(str(x) for x in spaced_kmer)
            seqs_kmers[s] = seqs_kmers.get(s, 0) + 1
            
            if i == 0 or i == len(sequence)-f:
                #solid_kmer = sequence[i:i+f]
                addon = 0
                if i == 0:
                    starts[s] = starts.get(s, 0) + 1
                else: 
                    ends[s] = ends.get(s, 0) + 1

            if s not in first_appearance:
                first_appearance[s] = all_startposis[index]
    return seqs_kmers, first_appearance, starts, ends



def unik(files, target_directory, has_reads=False, kmer_profile="1111110110110101110101011101011111101111"):
    profile = np.array([int(character) for character in kmer_profile])

    all_reads = []
    all_read_arr = []
    all_startposis = []
    all_startposis_arr = []
    if has_reads:
        for filename in files:
            reads, _ = misc.read_fasta(filename)
            all_reads += reads
            all_read_arr.append(reads)

    else:
        full_seqs = []
        for filename in files:
            reads, _ = misc.read_fasta(filename)
            full_seqs.append(reads)

        for seq in full_seqs:
            reads, startposis = misc.get_reads_from(seq, read_len=250, cov=20)
            all_reads = all_reads + reads
            all_read_arr.append(reads)
            all_startposis = all_startposis + startposis
            all_startposis_arr.append(startposis)


    all_seqs = []
    for i in range(len(all_reads)):
        seq = misc.parse_nucleotides(all_reads[i])
        all_seqs.append(seq)

    seqs_kmers, first_appearance, starts, ends = build_hashmaps(all_seqs, all_startposis, profile)
    all_raw, all_cc = cc.compute_counts(all_seqs, starts, ends, seqs_kmers, profile)
    uk = unique_kmers.compute_unique_kmers(all_seqs, all_raw, all_cc, profile)
    clusterIDs, counts, _ = clust.cluster(uk, len(all_seqs), first_appearance)
    misc.save_clusters(clusterIDs, all_reads, target_directory)
    