import utility 
import numpy as np
from importlib import reload
reload(utility)
import sys
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from multiprocessing import Pool

HOME = '/cluster/home/ajoudaki/projects2019-string-embedding/python_tools'

file_name = 'seqs0'
data_path = HOME + '/data/' + file_name + '.npz'
result_path = HOME + '/results/' + file_name + '.ann'
result_temp = HOME + '/results/temp_index'

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

proj_dim = 120
k_small = 2 
k_big = 100
num_trees = 3
num_neighbours = 40
num_neighbours_internal = 30
step_index = 1
step_search = step_index


kmers, s_kmer_vals, kmer_pos, kmer_seq_id = utility.list_kmers_simple(seqs, vals = vals,  
                                         k = k_small, addy = True, padding = int(k_big))
kmer_vals = utility.get_kmver_vals(s_kmer_vals, k_big)
convolved = utility.random_projection(kmers, k_big, proj_dim)
sids = np.array(kmer_seq_id)

def compute_index_search(si):
    index_indices = np.nonzero(sids==si)[0]
    search_indices = np.nonzero(sids!=si)[0]
    knn_index = utility.build_index(convolved, index_indices, num_trees)
    #NN,NN_dist = utility.knn_search_value(knn_index, convolved[search_indices,:],index_indices,  num_neighbours)
    NN,NN_dist = utility.knn_search_value(knn_index, convolved,search_indices,index_indices,  num_neighbours)
    np.savez(result_temp+str(si)+'.npz', NN=NN, NN_dist=NN_dist, search_indices = search_indices)

number_seqs = len(seqs)
with Pool() as pool:
    pool.map(compute_index_search, range(number_seqs))
    pool.close()
    pool.join()


# merge all NN indices
total_len = len(kmers)
NN = [[] for i in range(total_len)]
NN_dist = [[] for i in range(total_len)]
for sid in range(len(seqs)):
    Res = np.load(result_temp+str(sid)+'.npz')
    nn = Res['NN']
    nn_dist = Res['NN_dist']
    search_indices = Res['search_indices']
    for i,si in enumerate(search_indices):
        NN[si].extend(nn[i])
        NN_dist[si].extend(nn_dist[i])
#print('NN and NN_dist shape , len, len(NN[0]), len(NN[1]) ' )
#print(len(NN), len(NN[0]), len(NN[1]))
#print(len(NN_dist), len(NN_dist[0]), len(NN_dist[1]))
#for si in range(total_len):

search_indices = range(0,len(kmers),step_search)
quad = utility.eval_results(search_indices, NN, kmers, kmer_vals, kmer_pos, kmer_seq_id, num_seqs, options)

print('seq len = ', len(seqs[0]), ', num seqs = ', len(seqs))

print('hit ratio: ', len(quad)/len(search_indices)/num_neighbours)
