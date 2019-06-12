from collections import Counter
import numpy as np
import select
import sys
import argparse
from tqdm import tqdm
import time
from annoy import AnnoyIndex


def get_kmver_vals(s_kmer_vals, k_big):
    gene_indices = np.nonzero(np.convolve((np.array(s_kmer_vals)>0)*1,np.ones(k_big), mode='same')>0.5)[0]
    kmer_vals = np.zeros(shape=len(s_kmer_vals))
    for i in gene_indices:
        val,freq = Counter(s_kmer_vals[i-int(k_big/2):i+int(k_big/2)]).most_common(1)[0]
        kmer_vals[i] = val
    return kmer_vals

def my_print(string):
    print(string)
    sys.stdout.flush()
    

def eval_results(search_indices, NN, kmers, kmer_vals, kmer_pos, kmer_seq_id, num_seqs, options):
    t0 = time.time()
    quadrupples = []
    total_count = 0
    total_correct = 0
    
    total_len = len(search_indices)
    for kii in tqdm(range(total_len)):
        ki = search_indices[kii]
        kmer  = kmers[ki]
        kval = kmer_vals[ki]
        pos = kmer_pos[ki]
        sid = kmer_seq_id[ki]
        neighbour_indices = NN[kii]
        for ki2 in neighbour_indices:
            kmer2  = kmers[ki2]
            kval2 = kmer_vals[ki2]
            pos2 = kmer_pos[ki2]
            sid2 = kmer_seq_id[ki2]
            if sid!=sid2:
                total_count = total_count + 1
                if kval==kval2 and kval>0:
                    total_correct = total_correct + 1
                    s1,s2 = np.sort([sid2,sid])
                    quadrupples.append(((s1,s2),(kval, kval2)))

    quadrupples_set = set(quadrupples) 
#     quadrupples = sorted(quadrupples)
#     [print(uquad[i],', ', cquad[i]) for i in range(len(uquad))]
    Op = (options[()])
    print('evaluation time = ', time.time()-t0) 
    print('folls possitive = ', (total_count - total_correct)/total_count)
    Total = 0
    for num_seq in num_seqs:
        Total = Total + Op.num_seq*(Op.num_seq-1)/2 * Op.num_genes
    print('recall percentage ', len(quadrupples_set)/Total)
    return quadrupples


def all_kmers(Sigma, klen):
    unique_kmers = []
    sigma = len(Sigma)
    arr = [0]*(klen+1)
    while arr[klen]==0:
        s = ''.join([Sigma[arr[i]] for i in range(klen)])
        unique_kmers.append(s)
        j = 0
        arr[0] = arr[0] + 1
        while arr[j]>=sigma:
            arr[j] = 0
            j = j + 1
            arr[j] = arr[j] + 1
    return unique_kmers


def random_projection(skmers, k_big, proj_dim):
    np.random.seed(0)
    total_slen = len(skmers)
    #unique_skmers = list(set(skmers))
    unique_skmers = all_kmers('acgtxy',len(skmers[0]))
    nkmers = len(unique_skmers)

    kmer2id = {}
    for ki, kmer in enumerate(unique_skmers):
        kmer2id[kmer] = ki
    random_matrix = np.random.randint(0,2, size=(nkmers,proj_dim))*2 - 1

    my_print(' computing the projection')

    projection  = np.zeros(shape=(total_slen, proj_dim), dtype = np.int8)
    for ki in tqdm(range(total_slen)): 
        kid = kmer2id[skmers[ki]]
        projection[ki][:] = random_matrix[kid,:]

    my_print(' computing the convolution')

    convolved = np.zeros(shape = projection.shape, dtype = np.int8)
    for di in tqdm(range(proj_dim)):
        convolved[:,di] = np.convolve(projection[:,di],np.ones(k_big), mode='same')

    return convolved

def build_index(matrix, indices, num_trees, metric, index_path, verbose = True):
    total_len = len(indices)
    proj_dim = matrix.shape[1]
    # compute neighbors using annoy
    t0 = time.time()

    index = AnnoyIndex(proj_dim, metric= metric)  # Length of item vector that will be indexed
    index.on_disk_build(index_path) 
    for i in range(total_len):
        index.add_item(i, matrix[indices[i],:])
    index.build(num_trees)
    
    if verbose:
        my_print('time to build '+str(num_trees)+' trees = '+str(time.time()-t0)) 
    
    return index


def knn_search(index, indices, num_neighbours, verbose = True):
    if verbose:
        my_print('nearest neighbours search:') 
    t0 = time.time()
    
    if isinstance(indices[0], bool):
        indices = np.nonzero(indices)
    total_len = len(indices)

    NN = []
    NN_dist = []
    for kii in tqdm(range(total_len)):
        ki = indices[kii]
        nn,dist = index.get_nns_by_item(ki, num_neighbours, include_distances=True)
        NN.append(nn)
        NN_dist.append(dist)
    NN = np.array(NN)
    NN_dist = np.array(NN_dist)
        
    if verbose:
        my_print('time search = '+str(time.time()-t0)) 
    
    return NN, NN_dist

def knn_search_value(index, matrix, search_indices, build_indices, num_neighbours, verbose = True):
    if verbose:
        my_print('nearest neighbours search:') 
    t0 = time.time()
    
    total_len = len(search_indices)

    NN = np.zeros(shape=(total_len,num_neighbours), dtype = np.int64)
    NN_dist = np.zeros(shape=(total_len,num_neighbours), dtype=np.float64)
    for ki in tqdm(range(total_len)):
        si = search_indices[ki]
        nn,nn_dist = index.get_nns_by_vector(matrix[si,:], num_neighbours, include_distances=True)
        NN[ki,:] = build_indices[nn]
        NN_dist[ki,:] = nn_dist
        
    if verbose:
        my_print('time search = '+str(time.time()-t0)) 
    
    return NN, NN_dist


def list_kmers(seqs, vals = [],  k = 10,  to_sort = False, addx = True, addy = False, unique = False,  print_kmers = False):
    kmers = []
    kmer_vals = []
    

    for i, s in enumerate(seqs):
        if addx:
            if not addy:
                S = ('x'*k + s + 'x'*k)
            else:
                S = ('x'*k + s + 'y'*k)
            if len(vals)!=0:
                Val = [0]*k
                Val.extend(vals[i])
                Val.extend([0]*k)

        else: 
            S = s 
            Val = vals[i]

        for i in range(len(S)-k+1):
            kmers.append(S[i:i+k])
            val,freq = Counter(Val[i:i+k]).most_common(1)[0]

            if val>0 and freq>k*.7:
                kmer_vals.append(val)
            else:
                kmer_vals.append(0)

    if unique:
        kmers = list(set(kmers))
    Map = {}
    for kmer,kmer_val in zip(kmers, kmer_vals):
        if kmer not in Map.keys():
            Map[kmer] = [kmer_val]
        else:
            Map[kmer].append(kmer_val)
            kmers = []
            kmer_vals = []

        for k,v in Map.items():
            kmers.append(k)
            kmer_vals.append(Counter(v).most_common(1)[0][0])


    if print_kmers:
        for s in kmers:
            print(s)

    if to_sort:
        kmers = sorted(kmers)


    return kmers, kmer_vals

def list_kmers_simple(seqs, vals ,  k = 10,  to_sort = False, 
                      addx = True, addy = False, 
                      unique = False,  print_kmers = False, 
                      threshold = .7,
                      simple_val = True,
                      padding = -1):
    
    my_print('computing the list of kmers and values ')
    
    kmers = []
    kmer_vals = []
    kmer_seq_id = []
    kmer_pos = []
    
    if padding == -1:
        padding = k

    for si in tqdm(range(len(seqs))):
        s = seqs[si]
        if addx:
            if not addy:
                S = ('x'*padding + s + 'x'*padding)
            else:
                S = ('x'*padding + s + 'y'*padding)
            if len(vals)!=0:
                Val = [0]*padding
                Val.extend(vals[si])
                Val.extend([0]*padding)

        else: 
            S = s 
            Val = vals[si]
        for i in range(len(S)-k+1):
            kmers.append(S[i:i+k])
            if simple_val:
                kmer_vals.append(Val[i+int(k/2)])
            else:
                val,freq = Counter(Val[i:i+k]).most_common(1)[0]

                if val>0 and freq>k*(threshold):
                    kmer_vals.append(val)
                else:
                    kmer_vals.append(0)
        kmer_seq_id.extend([si]*(len(S)-k+1))
        kmer_pos.extend([i for i in range(len(S)-k+1)])
                
    return kmers, kmer_vals, kmer_pos, kmer_seq_id



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kmer",type=int, default = 10, help="small kmer size to be counted a")
    parser.add_argument("-x", "--addx", default = False, action="store_true", help="small kmer size to be counted a")
    parser.add_argument("-y", "--addy", default = False, action="store_true", help="small kmer size to be counted a")
    parser.add_argument("-s", "--sort", default = False, action="store_true", help="small kmer size to be counted a")
    parser.add_argument("-u", "--unique", default = False, action="store_true", help="small kmer size to be counted a")
    args = parser.parse_args()

    seqs = []

    for line in sys.stdin:
        seqs.append(line[:-1])

    list_kmers(seqs, args.kmer, args.sort, args.addx, args.addy, args.unique,  print_kmers= True)
