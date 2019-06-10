import utility 
import numpy as np
from importlib import reload
reload(utility)
import sys, os, shutil
import time
from tqdm import tqdm
import matplotlib.pyplot as pl
from optparse import OptionParser
from annoy import AnnoyIndex

#HOME = '/cluster/home/ajoudaki/projects2019-string-embedding/python_tools'
HOME = '/cluster/work/grlab/share/databases/genomes/synthetic'

def load_paths(file_name):
    data_path = HOME + '/data/' + file_name + '.npz'
    result_path = HOME + '/results/' + file_name 
    index_path = result_path + '/index/'
    search_path = result_path + '/search/'
    eval_path = result_path + '/eval/'
    npz_path = result_path + '/npz/'

    return result_path, index_path, search_path, eval_path, npz_path

def load_files(file_name, clean):
    data_path = HOME + '/data/' + file_name + '.npz'
    result_path = HOME + '/results/' + file_name 
    index_path = result_path + '/index/'
    search_path = result_path + '/search/'
    eval_path = result_path + '/eval/'
    npz_path = result_path + '/npz/'

    paths = (index_path, search_path, npz_path, eval_path)

    if clean==True:
        if os.path.exists(result_path):
            shutil.rmtree(result_path)
        os.makedirs(result_path)
        return

    for p in paths:
        if not os.path.exists(p):
            os.makedirs(p)

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
    print('num seqs = ', len(seqs), ' mean lenth of seqs = ', np.mean([len(seqs[i]) for i in range(len(seqs))]))
    print(Op)
    
    return seqs, vals, num_seqs, options, Op, paths
    


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f','--file-name',dest='file_name', default = 'seqs0', type='string', help = 'seq file to load')
    parser.add_option('-i','--seq-id',dest='seq_id', default = 0, type='int', help = 'seq id to build index for')
    parser.add_option('-j','--search-seq-id',dest='search_seq_id', default = 0, type='int', help = 'search sequence id')
    parser.add_option('-k','--small-k',dest='k_small', default = 3, type='int', help = 'small k to do kmer counting')
    parser.add_option('-K','--big-k',dest='k_big', default = 100, type='int', help = 'big K to slide over the strings')
    parser.add_option('-D','--proj-dim',dest='proj_dim', default = 120, type='int', help = 'number of dimensions to project onto')
    parser.add_option('-n','--num-trees',dest='num_trees', default = 2, type='int', help = 'number of trees to be built')
    parser.add_option('-s','--step-build',dest='step_build', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-m','--metric',dest='metric', default = 'euclidean', type='string', help = 'metric to be used')
    parser.add_option('-S','--step-search',dest='step_search', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-N','--num-neighbors',dest='num_neighbors', default = 50, type='int', help = 'number of nearest neighbors to search')
    parser.add_option('-M','--memory',dest='memory', default = 30, type='int', help = 'memory to allocate for each process')
    parser.add_option('-t','--target',dest='target', default = 'build', type='string', help = 'metric to be used')
    parser.add_option('-T','--forward-target',dest='forward_target', default = 'all', type='string', help = 'metric to be used')
    (options, args) = parser.parse_args()

    sid = options.seq_id
    sid2 = options.search_seq_id
    proj_dim = options.proj_dim
    k_small = options.k_small 
    k_big = options.k_big
    num_trees = options.num_trees
    step_build = options.step_build
    step_search = options.step_search
    num_neighbors = options.num_neighbors
    metric = options.metric


    if option.target=='all':
        pass


    elif options.target=='clean':
        load_files(options.file_name, clean=True)

    elif options.target=='build':
    
        seqs, vals, num_seqs, opts, Op, paths = load_files(options.file_name, clean=False)
        index_path, search_path, npz_path , eval_path= paths
        
        kmers, s_kmer_vals, kmer_pos, kmer_seq_id = utility.list_kmers_simple([seqs[sid]], vals = [vals[sid]],  
                                                 k = k_small, addy = True, padding = int(k_big))
        kmer_vals = utility.get_kmver_vals(s_kmer_vals, k_big)
        build_indices = np.arange(0,len(kmers),step_build)

        convolved = utility.random_projection(kmers, k_big, proj_dim)
        np.savez(npz_path+'data'+str(sid)+'.npz', 
                 convolved=convolved, 
                 num_seqs = num_seqs,
                 kmers = kmers, 
                 kmer_pos = kmer_pos,
                 kmer_vals = kmer_vals, 
                 build_indices = build_indices)

        knn_index = utility.build_index(matrix = convolved, 
                                        indices = build_indices, 
                                        num_trees = num_trees, 
                                        metric = metric, 
                                        index_path = index_path + str(sid) + '.ann')
        np.savez(index_path+str(sid)+'.done', sid=sid)

    elif options.target=='search':
        #seqs, vals, num_seqs, opts, Op, paths = load_files(options.file_name, clean=False)
        #index_path, search_path, npz_path , eval_path= paths
        result_path, index_path, search_path, eval_path, npz_path = load_paths(options.file_name)

        data1 = np.load(npz_path+'data'+str(sid)+'.npz')
        knn_index = AnnoyIndex(proj_dim, metric)
        knn_index.load(index_path + str(sid) + '.ann', prefault=True)
        build_indices = data1['build_indices']

        data2 = np.load(npz_path+'data'+str(sid2)+'.npz')
        convolved = data2['convolved']
        search_indices = np.arange(0,convolved.shape[0],step_search)
        NN, NN_dist = utility.knn_search_value(knn_index, convolved, search_indices, build_indices,num_neighbors)
        np.savez(search_path + str(sid)+'_'+str(sid2) + '.npz', NN=NN, NN_dist=NN_dist, search_indices = search_indices)
        #np.savez(search_path+str(sid)+'_'+str(sid2)+'.done', sid=sid)
        print('done',file=open(search_path +str(sid)+'_'+str(sid2)+'.done','w+'))

    elif options.target=='eval':
        seqs, vals, num_seqs, opts, Op, paths = load_files(options.file_name, clean=False)
        index_path, search_path, npz_path , eval_path= paths

        data1 = np.load(npz_path+'data'+str(sid)+'.npz')
        data2 = np.load(npz_path+'data'+str(sid2)+'.npz')
        kvals1 = data1['kmer_vals']
        kvals2 = data2['kmer_vals']
        search_data = np.load(search_path + str(sid)+'_'+str(sid2) + '.npz')
        NN = search_data['NN']
        NN_dist = search_data['NN_dist']
        search_indices = search_data['search_indices']
        build_indices = data1['build_indices']

        total_len = len(search_indices)
        total_count = 0
        total_correct = 0
        quadrupples = []
        for kii in tqdm(range(total_len)):
            ki = search_indices[kii]
            kval2 = kvals2[ki]
            neighbor_indices = NN[kii]
            for ki2 in neighbor_indices:
                kval1 = kvals1[ki2]
                total_count = total_count + 1
                if kval1==kval2 and kval1>0:
                    total_correct = total_correct + 1
                    s1,s2 = np.sort([sid,sid2])
                    quadrupples.append( ( (s1,s2) , (kval1, kval2) ) )
        quadrupples = set(quadrupples)
        np.savez(eval_path +str(sid)+'_'+str(sid2)+'.npz', quadrupples=quadrupples)
        #np.savez(eval_path +str(sid)+'_'+str(sid2)+'.done', sid=sid)
        print('done',file=open(eval_path +str(sid)+'_'+str(sid2)+'.done','w+'))
        print('false positive = ', (total_count - total_correct)/total_count)
        print('recall = ', len(quadrupples)*1.0/Op.num_genes)

