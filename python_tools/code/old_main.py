import utility 
import numpy as np
from importlib import reload
reload(utility)
import sys
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from multiprocessing import Pool
from annoy import AnnoyIndex


def load_data(seq_number):
    file_name = 'seqs' + str(seq_number)
    data_path = HOME + '/data/' + file_name + '.npz'
    Res = np.load(data_path)
    seqs = Res['seqs']
    vals = Res['vals']
    options = Res['options']
    Op = (options[()])
    if int(file_name[4:])>=7:
        gene_lens = Res['gene_lens']
        num_seqs = Res['num_seqs']
    else:
        num_seqs = [Op.num_seq]*Op.repeat
        
    print(Op)
    print('seq len = ', len(seqs[0]))
    print('num seqs = ', len(seqs))
    return seqs, vals, options, num_seqs

def compute_kmers(seqs, vals, k_small, k_big):
    kmers, s_kmer_vals, kmer_pos, kmer_seq_id = utility.list_kmers_simple(seqs, vals = vals,  
                                             k = k_small, addy = True, padding = int(k_big))
    kmer_vals = utility.get_kmver_vals(s_kmer_vals, k_big)
    kmer_sids = np.array(kmer_seq_id)
    return kmers, kmer_vals, kmer_pos, kmer_sids


def compute_index(si, kmer_sids, num_trees, convolved, result_path):
    build_indices = np.nonzero(kmer_sids==si)[0]
    knn_index = utility.build_index(convolved, build_indices, num_trees)
    knn_index.save(result_path + 'index_' + str(si) + '.ann')

def search_index(si, kmer_sids, convolved, num_neighbours, chunk_i, chunk_num, result_path):
    project_dim = convolved.shape[1]
    knn_index = AnnoyIndex(project_dim, metric = "euclidean")
    knn_index.load(result_path + 'index_' + str(si) + '.ann', prefault = True)
    nn, nn_dist = knn_index.get_nns_by_vector(convolved[0], 10, include_distances = True)
    build_indices = np.nonzero(kmer_sids==si)[0]
    search_indices = np.nonzero(kmer_sids!=si)[0]
    chunk_size = len(search_indices)*1.0/chunk_num
    chunk_size = int(np.ceil(chunk_size))
    search_indices = search_indices[chunk_i*chunk_size:(chunk_i+1)*chunk_size]
    NN,NN_dist = utility.knn_search_value(knn_index, convolved,search_indices,build_indices,  num_neighbours)
    np.savez(result_path+ 'search_s'+str(si)+'_c'+ str(chunk_i)+'.npz', NN=NN, NN_dist=NN_dist, search_indices = search_indices)



if __name__ == '__main__':
    HOME = '/cluster/home/ajoudaki/projects2019-string-embedding/python_tools'

    seq_number = 0
    file_name = 'seqs' + str(seq_number)
    result_path = HOME + '/results/' + file_name + '/'
    #result_temp = HOME + '/results/temp_result'

    proj_dim = 120
    k_small = 2 
    k_big = 100
    num_trees = 1
    num_neighbours = 40
    num_neighbours_internal = 30
    chunk_num = 2
    step_index = 1
    step_search = step_index


    seqs, vals, options, num_seqs = load_data(seq_number)
    kmers, kmer_vals, kmer_pos, kmer_sids = compute_kmers(seqs, vals, k_small, k_big)
    convolved = utility.random_projection(kmers, k_big, proj_dim)


    number_seqs = len(seqs)
    for si in range(number_seqs):
        compute_index(si, kmer_sids, num_trees, convolved, 
                result_path = result_path)
    for si in range(number_seqs):
        for ci in range(chunk_num):
            search_index(si, kmer_sids, convolved, 
                    num_neighbours, 
                    chunk_i = ci, chunk_num = chunk_num, 
                    result_path = result_path)

    # merge all NN indices
    total_len = len(kmers)
    merge_chunk_num = 10
    merge_chunk_size = total_len*1.0/merge_chunk_num
    merge_chunk_size = int(np.ceil(merge_chunk_size))


    # merge in chunks of the total size
    NN = [[] for i in range(total_len)]
    NN_dist = [[] for i in range(total_len)]
    for sid in range(len(seqs)):
        for ci in range(chunk_num):
            Res = np.load(result_path+ 'search_s'+str(sid)+'_c'+str(ci)+'.npz')
            nn = Res['NN']
            nn_dist = Res['NN_dist']
            search_indices = Res['search_indices']
            for i,si in enumerate(search_indices):
                NN[si].extend(nn[i])
                NN_dist[si].extend(nn_dist[i])

    search_indices = range(0,len(kmers),step_search)
    quad = utility.eval_results(search_indices, NN, kmers, kmer_vals, kmer_pos, kmer_sids, num_seqs, options)

    print('seq len = ', len(seqs[0]), ', num seqs = ', len(seqs))

    print('hit ratio: ', len(quad)/len(search_indices)/num_neighbours)

