import utility 
import numpy as np
from importlib import reload
reload(utility)
import sys, os, time
from tqdm import tqdm
import matplotlib.pyplot as plt

HOME = '/cluster/home/ajoudaki/projects2019-string-embedding/python_tools'

file_name = 'seqs0'
data_path = HOME + '/data/' + file_name + '.npz'
results_path = HOME + '/results/' 
temp_dir = results_path + file_name 
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

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

proj_dim = 120
k_small = 2 
k_big = 100
num_trees = 3
num_neighbours = 200


kmers, s_kmer_vals, kmer_pos, kmer_seq_id = utility.list_kmers_simple(seqs, vals = vals,  
                                         k = k_small, addy = True, padding = int(k_big))
kmer_vals = utility.get_kmver_vals(s_kmer_vals, k_big)
search_indices = np.array([i for i in range(0,len(kmers),1)])
build_indices = np.array([i for i in range(0,len(kmers),1)])

convolved = utility.random_projection(kmers, k_big, proj_dim)

knn_index = utility.build_index(convolved, build_indices, num_trees)
NN,NN_dist = utility.knn_search_value(knn_index, convolved, search_indices, build_indices, num_neighbours)
    
quad = utility.eval_results(search_indices, NN, kmers, kmer_vals, kmer_pos, kmer_seq_id, num_seqs, options)

print('seq len = ', len(seqs[0]), ', num seqs = ', len(seqs))

print('hit ratio: ', len(quad)/len(search_indices)/num_neighbours)
