import utility 
import numpy as np
from importlib import reload
reload(utility)
import sys, os, shutil, time
from tqdm import tqdm
import matplotlib.pyplot as pl
from optparse import OptionParser
from annoy import AnnoyIndex

HOME = '/cluster/work/grlab/share/databases/genomes/synthetic'

def load_paths(file_name):
    data_path = HOME + '/data/' + file_name + '.npz'
    result_path = HOME + '/results/' + file_name +'/'  
    index_path = result_path + 'index/'
    search_path = result_path + 'search/'
    eval_path = result_path + 'eval/'
    npz_path = result_path + 'npz/'
    log_path = result_path + 'log/'

    return data_path, result_path, index_path, search_path, eval_path, npz_path, log_path

def load_files(file_name, clean):
    data_path, result_path, index_path, search_path, eval_path, npz_path, log_path = load_paths(file_name)
    paths = (index_path, search_path, npz_path, eval_path, log_path)

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
    
    return seqs, vals, num_seqs, options, Op
    

def get_job_path(job, job_paths):
    S = job_paths[job[0]] + job[0]
    for i in range(1,len(job)):
        S = S + '_' + str(job[i])
    return S + '.done'

def get_job_cmd(job):
    S = ' -t ' + job[0] 
    argnames = ['', ' -i ', ' -j '] 
    for i in range(1,len(job)):
        S = S + argnames[i] + str(job[i])
    return S 

def add_job(Map, job, deps):
    Map[job] = deps

def crawl_dep_recurse(Map, Set, curr):
    if curr not in Set:
        Set.add(curr)
        if curr in Map.keys():
            for job in Map[curr]:
                crawl_dep_recurse(Map, Set, job)

def crawl_dep(Map, curr):
    Set = set()
    crawl_dep_recurse(Map, Set,curr)
    #return list(Set)
    return Set


def run_jobs(jobs_flat, job_dep, command, job_paths):
    _, _, _, _, _, _, log_path = load_paths(file_name)
    started = dict()
    for j in jobs_flat:
        started[j] = False
    fcommands = open(log_path + 'commands.sh', 'w+')
    final_job = False
    while not final_job:
        time.sleep(1)
        for job in jobs_flat:
            jp = get_job_path(job,job_paths)
            if not os.path.exists(jp) and not started[job]:
                deps_satisfied = [os.path.exists(get_job_path(j,job_paths)) for j in job_dep[job] ] 
                if all(deps_satisfied):
                    cmd_combo = command + get_job_cmd(job)  
                    if job[0]=='final':
                        final_job = True
                        os.system('mv lsf* '+ log_path)
                        with open(log_path + 'final.result', 'r') as results:
                            print(results.read())
                    else:
                        os.system(cmd_combo)
                    started[job] = True
                    print(cmd_combo)
                    print(cmd_combo, file = fcommands)
    fcommands.close()



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f','--file-name',dest='file_name', default = 'seqs0', type='string', help = 'seq file to load')
    parser.add_option('-i','--seq-id',dest='seq_id', default = 0, type='int', help = 'seq id to build index for')
    parser.add_option('-j','--search-seq-id',dest='search_seq_id', default = 0, type='int', help = 'search sequence id')
    parser.add_option('-k','--small-k',dest='small_k', default = 3, type='int', help = 'small k to do kmer counting')
    parser.add_option('-K','--big-k',dest='big_k', default = 100, type='int', help = 'big K to slide over the strings')
    parser.add_option('-D','--proj-dim',dest='proj_dim', default = 120, type='int', help = 'number of dimensions to project onto')
    parser.add_option('-n','--num-trees',dest='num_trees', default = 2, type='int', help = 'number of trees to be built')
    parser.add_option('-s','--step-build',dest='step_build', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-m','--metric',dest='metric', default = 'euclidean', type='string', help = 'metric to be used')
    parser.add_option('-S','--step-search',dest='step_search', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-N','--num-neighbors',dest='num_neighbors', default = 50, type='int', help = 'number of nearest neighbors to search')
    parser.add_option('-M','--memory',dest='memory', default = 30, type='int', help = 'memory to allocate for each process')
    parser.add_option('-t','--target',dest='target', default = 'build', type='string', help = 'metric to be used')
    parser.add_option('-T','--forward-target',dest='forward_target', default = 'merge', type='string', help = 'metric to be used')
    (options, args) = parser.parse_args()

    sid = options.seq_id
    sid2 = options.search_seq_id
    proj_dim = options.proj_dim
    k_small = options.small_k 
    k_big = options.big_k
    num_trees = options.num_trees
    step_build = options.step_build
    step_search = options.step_search
    num_neighbors = options.num_neighbors
    metric = options.metric
    file_name = options.file_name

    data_path, result_path, index_path, search_path, eval_path, npz_path, log_path = load_paths(file_name)
    job_paths = {'clean':result_path, 
            'build': index_path,
            'search' : search_path,
            'eval': eval_path,
            'merge': result_path ,
            'final': result_path }


    if options.target=='all':
        summary = np.load(result_path+'num.npz')
        N = summary['N']
        jobs_dep = dict()
        add_job(jobs_dep,('clean',),[])
        build_jobs = []
        for i in range(N):
            job = ('build',i) 
            add_job(jobs_dep,job,[('clean',)])
            build_jobs.append(job)
        search_jobs = []
        for i in range(N):
            for j in range(N):
                if i!=j:
                    job = ('search',i,j)
                    add_job( jobs_dep,job,[('build',i),('build',j)]  )
                    search_jobs.append(search_jobs)
        eval_jobs = []
        for i in range(N):
            for j in range(N):
                if i!=j:
                    job = ('eval',i,j)
                    add_job( jobs_dep, job, [ ('search',i,j) ]  )
                    eval_jobs.append(job)
        add_job(jobs_dep, ('merge',), eval_jobs)

        dep_job = {'final':[], 
                'clean': [('clean',)],
                'build' : build_jobs,
                'search': search_jobs,
                'eval': eval_jobs, 
                'merge': [('merge',)] }
        add_job(jobs_dep,('final',), dep_job[options.forward_target]) 
        jobs_flat = crawl_dep(jobs_dep, ('final',))
        command = 'bsub -R "rusage[mem=' + str(options.memory) + '000]"  python bk_worker_main.py '

        for k,v in options.__dict__.items():
            if k not in ['target', 'seq_id', 'search_seq_id', 'target', 'forward_target']:
                command = command + ' --' + k.replace('_','-') + ' ' + str(v) 
        run_jobs(jobs_flat, jobs_dep, command, job_paths) 

    elif options.target=='final':
        with open(log_path + 'final.result', 'r') as results:
            print(results.read())
        os.system('mv lsf* '+ log_path)
        job_done = ('final', )

        
    elif options.target=='merge':
        quads = set()
        summary = np.load(result_path+'num.npz')
        N = summary['N']
        num_seqs = summary['num_seqs']
        Op = summary['Op']
        Op = Op[()]
        for i in range(N):
            for j in range(N):
                if i!=j:
                    res = np.load(eval_path +str(i)+'_'+str(j)+'.npz') 
                    q = res['quadrupples']
                    quads = quads.union(q[()])
        Total = 0
        for ns in num_seqs:
            Total = Total + ns*(ns-1)/2 * Op.num_genes
        fresult = open(log_path + 'final.result', 'w+')
        print('#'*50, file=fresult)
        print('final recall : ', len(quads)*1.0/Total, file=fresult)
        print('#'*50, file = fresult)
        fresult.close()

        job_done = ('merge', )

    elif options.target=='clean':
        load_files(options.file_name, clean=True)
        seqs, vals, num_seqs, opts, Op = load_files(file_name, clean=False)
        np.savez(result_path+'num.npz', N=len(seqs), num_seqs = num_seqs, options = options, Op = Op)

        job_done = ('clean', )

    elif options.target=='build':
    
        seqs, vals, num_seqs, opts, Op  = load_files(file_name, clean=False)
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
        job_done = ('build', sid)

    elif options.target=='search':

        data1 = np.load(npz_path+'data'+str(sid)+'.npz')
        knn_index = AnnoyIndex(proj_dim, metric)
        knn_index.load(index_path + str(sid) + '.ann', prefault=True)
        build_indices = data1['build_indices']

        data2 = np.load(npz_path+'data'+str(sid2)+'.npz')
        convolved = data2['convolved']
        search_indices = np.arange(0,convolved.shape[0],step_search)
        NN, NN_dist = utility.knn_search_value(knn_index, convolved, search_indices, build_indices,num_neighbors)
        np.savez(search_path + str(sid)+'_'+str(sid2) + '.npz', NN=NN, NN_dist=NN_dist, search_indices = search_indices)
        job_done = ('search', sid, sid2)

    elif options.target=='eval':

        data1 = np.load(npz_path+'data'+str(sid)+'.npz')
        data2 = np.load(npz_path+'data'+str(sid2)+'.npz')
        summary = np.load(result_path + 'num.npz')
        Op = summary['Op']
        Op = Op[()]
        num_seqs = data1['num_seqs']
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
        np.savez(eval_path +str(sid)+'_'+str(sid2)+'.npz', 
                quadrupples = quadrupples, 
                total_count = total_count, 
                total_correct = total_correct)

        print('false positive = ', (total_count - total_correct)/total_count)
        print('recall = ', len(quadrupples)*1.0/Op.num_genes)
        
        job_done = ('eval', sid, sid2)

    # print done 
    if options.target!='all':
        Job_path = get_job_path(job_done, job_paths)
        print('dummy', file=open(Job_path, 'w+'))

