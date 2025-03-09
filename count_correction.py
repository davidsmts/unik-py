import misc


def count_correction(target, starts, ends, hashmap, profile):
    f = len(profile)
    counts = []
    corrected_counts = []
    raw_counts = []
    maximum = 0
    mci = 0

    # create all sk-mers
    skmers = []
    for i in range(len(target)-f+1):
        skmer = misc.get_skmer(target[i:i+f], profile)
        skmers.append(skmer)
        
    for i in range(len(target)-f+1):
        K_i = skmers[i]
        c_K_i = hashmap[K_i]
        counts.append(c_K_i)
        
        startcount = 0
        endcount = 0
        for j in range(0, i+1):
            K_j = skmers[j]
            #c_K_j = hashmap[K_j]
            #if j != max(i-f+1, 0):
            startcount += starts.get(K_j, 0)
            if j != i:
                endcount += ends.get(K_j, 0)
        k_dash_i = c_K_i - startcount + endcount
        corrected_counts.append(k_dash_i)
    raw_counts = []
    
    return counts, corrected_counts


def compute_counts(all_seqs, starts, ends, seqs_kmers, profile):
    # COMPUTE MAXCOUNTS
    all_corr_counts = []
    all_raw = []
    for j in range(len(all_seqs)):
        if j%10000 == 0:
            print(j)
        target = all_seqs[j]
        raw, corr_counts = count_correction(target, starts, ends, seqs_kmers, profile)
        all_raw.append(raw)
        all_corr_counts.append(corr_counts)

    return all_raw, all_corr_counts


