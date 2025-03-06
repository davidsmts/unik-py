import misc
import numpy as np
import sys

def get_candidates(seq, unique_positions, profile):
    candidates = []
    candidate_pos = []
    for pos in unique_positions:
        candidates.append(misc.get_skmer(seq[pos:pos+len(profile)], profile))
        candidate_pos.append(pos)
    return candidates, candidate_pos


def get_unique_kmers_from_read(k, kcorr, L=250, f=40, smin=0, delta=1):
    uniq_cov = sys.maxsize
    nonuniq_cov = 0
    for i in range(1, L-f+1):
        if kcorr[i] - kcorr[i-1] >= smin and kcorr[i] < uniq_cov and smin <= k[i]:
            uniq_cov = k[i-1]
            nonuniq_cov = kcorr[i]
        if kcorr[i-1] - kcorr[i] >= smin and kcorr[i-1] < uniq_cov and smin <= k[i-1]:
            uniq_cov = k[i]
            nonuniq_cov = kcorr[i-1]

    uniques = []
    for i in range(L-f+1):
        if nonuniq_cov == 0 or (k[i] <= uniq_cov + np.sqrt(delta*uniq_cov) and kcorr[i] <= nonuniq_cov - smin and k[i] >= smin):
            uniques.append(i)

    return uniques


def compute_unique_kmers(all_seqs, all_raw, all_corr_counts, profile):
    #k, kcorr, L, f, smin, delta
    unique_kmers = []
    debug_unique_kmers = []
    cnt = 0
    m = 40
    for i in range(len(all_raw)):
        if i%10000 == 0:
            print(i)
        posis = get_unique_kmers_from_read(all_raw[i], all_corr_counts[i], L=250, f=40, smin=2, delta=2)
        if posis == []:
            print("Empty!")
            print(i)
            
        if posis != [j for j in range(len(all_raw[i])-39)]:
            cnt += 1
        seq = all_seqs[i]
        candidates, _ = get_candidates(seq, posis, profile)
        hashed = []
        for candidate in candidates:
            #tohash = ''.join([str(el) for el in candidate])
            tohash = ''.join([str(el) for el in candidate])
            hashed.append(mmh3.hash(tohash))
        sorted_indices = np.argsort(hashed)
        representatives = np.take_along_axis(np.array(candidates), sorted_indices, axis=0)
        for repres in representatives:
            unique_kmers.append([i, str(repres)])
        debug_unique_kmers.append([i, posis])
    #print(candidates)
    print(unique_kmers[:10])
    print(cnt)
    print(cnt / len(all_raw))
    return unique_kmers


