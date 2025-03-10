import numpy as np
import matplotlib.pyplot as plt
import mmh3
from collections import defaultdict, deque
import os.path as op
import misc
import unique_kmers
import count_correction as cc
import cluster as clust
import time
from pathlib import Path

def build_hashmaps(all_seqs, all_startposis, profile):
    # BUILD THE HASHMAP
    f = len(profile)
    starts = {}
    ends = {}
    first_appearance = {}
    seqs_kmers = {}
    print("Building hashmap")
    prev_index = 0
    print(len(all_seqs))
    for (index, sequence) in zip(range(len(all_seqs)), all_seqs):
        if index-prev_index > len(all_seqs)*0.1:
            prev_index = index
            print("hashmap progress: " + str(np.round(100*index/len(all_seqs))) + "%")
        for i in range(len(sequence) - f +1):
            skmer = misc.get_skmer(sequence[i:i+f], profile)
            seqs_kmers[skmer] = seqs_kmers.get(skmer, 0) + 1
            
            if i == 0 or i == len(sequence)-f:
                #solid_kmer = sequence[i:i+f]
                addon = 0
                if i == 0:
                    if skmer in starts:
                        starts[skmer] += 1
                    else:
                        starts[skmer] = 0
                else: 
                    if skmer in ends:
                        ends[skmer] += 1
                    else:
                        ends[skmer] = 0
                    #ends[skmer] = ends.get(skmer, 0) + 1
            if len(all_startposis) > 0:
                if skmer not in first_appearance:
                    first_appearance[skmer] = all_startposis[index]
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
            if filename != ".DS_Store":
                reads, _ = misc.read_fasta(filename)
                full_seqs += reads

        for seq in full_seqs:
            reads, startposis = misc.get_reads_from(seq, read_len=250, cov=7)
            all_reads = all_reads + reads
            all_read_arr.append(reads)
            all_startposis = all_startposis + startposis
            all_startposis_arr.append(startposis)

    all_seqs = []
    for i in range(len(all_reads)):
        seq = misc.parse_nucleotides(all_reads[i])
        all_seqs.append(seq)

    timestamp = time.time()
    tardir = "./tmp/" + str(timestamp) + "/"
    Path(tardir).mkdir(parents=True, exist_ok=True)
    seqs_kmers, first_appearance, starts, ends = build_hashmaps(all_seqs, all_startposis, profile)
    misc.save_dicts(tardir + "hashmaps", [seqs_kmers, first_appearance, starts, ends])
    all_raw, all_cc = cc.compute_counts(all_seqs, starts, ends, seqs_kmers, profile)
    #np.savetxt(tardir + "allraw.csv", np.array(all_raw), delimiter=",")
    #np.savetxt(tardir + "allcc.csv", np.array(all_cc), delimiter=",")
    uk = unique_kmers.compute_unique_kmers(all_seqs, all_raw, all_cc, profile)
    np.save(tardir + "uk.npy", uk)
    clusterIDs, counts, _ = clust.cluster(uk, len(all_seqs), first_appearance)
    misc.save_clusters(clusterIDs, all_reads, target_directory)
    