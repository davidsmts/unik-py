import numpy as np


def clusters_from_ids(clusters, all_reads):
    diff_clusts = 0
    clusts = {}
    for (i, el) in zip(range(len(clusters)), clusters):
        if el not in clusts:
            clusts[el] = [all_reads[i]]
            diff_clusts += 1
        else:
            clusts[el].append(i)

    clust_arr = []
    for key in clusts.keys():
        cluster = []
        for id in clusts[key]:
            cluster.append(all_reads[id])
        clust_arr.append(cluster)

    return clust_arr


def followToLocalMinimum(readIDs, clusterIDs):
    minIDs = []
    for i in range(len(readIDs)):
        id = readIDs[i]
        while id != clusterIDs[id]:
            id = clusterIDs[id]
        minIDs.append(id)
    return minIDs



def cluster(unique_kmers, num_reads, startposis):
    print("clustering")
    sortedarr = sorted(unique_kmers, key=lambda x: x[1])
    print("finished sort")
    clusterIDs = [i for i in range(num_reads)]
    reads_ = []
    #for i in range(len(clusterIDs)):
    #    clusterIDs[sortedarr[i][0]] = sortedarr[i][0]
    # identify index range and set to minimum read index
    # since the relative order of the elements is maintained, we will
    counts = [0 for _ in range(60)]
    idx0 = 0
    prev = 0
    while idx0 < len(sortedarr):
        #print("idx0 = " + str(idx0))
        curr_kmer = sortedarr[idx0][1]
        indices = [sortedarr[idx0][0]]
        idx1 = idx0 + 1
        if idx1 >= len(sortedarr):
            break
        while sortedarr[idx1][1] == curr_kmer:
            indices.append(sortedarr[idx1][0])
            # erste k-mer occurence nehmen
            if len(list(startposis.keys())) > 0:
                if np.abs(startposis[sortedarr[idx0][1]] - startposis[sortedarr[idx1][1]]) > 220:
                    reads_.append((sortedarr[idx0][0],sortedarr[idx1][0]))
            idx1 += 1
            if idx1 == len(sortedarr):
                break
        
        minIDs = followToLocalMinimum(indices, clusterIDs)
        minID = min(minIDs)
        for id in minIDs:
            clusterIDs[id] = minID
        for id in indices:
            clusterIDs[id] = minID
    
        idx0 = idx1
        prev = idx0
    
    #print(clusterIDs[:1000])
    clusterIDs = followToLocalMinimum([i for i in range(num_reads)], clusterIDs)
    return clusterIDs, counts, reads_