import numpy as np
import sys, os, shutil, time, subprocess, glob, re
from tqdm import tqdm
import matplotlib.pyplot as pl
from optparse import OptionParser
from annoy import AnnoyIndex
from attrdict import AttrDict
from importlib import reload
from datetime import datetime
import utility 
reload(utility)

# set the directory where the sequence data is located
PROJ_DIR = '/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic'
RESULT_DIR = '__UNASSIGNED__'

def load_paths(file_name):
    data_path = PROJ_DIR + '/data/' + file_name + '.npz'
    result_path = PROJ_DIR + '/results/' + RESULT_DIR +  '/'   
    index_path = result_path + 'index/'
    search_path = result_path + 'search/'
    eval_path = result_path + 'eval/'
    npz_path = result_path + 'npz/'
    log_path = result_path + 'log/'

    paths = (result_path, index_path, search_path, npz_path, eval_path, log_path)
    for p in paths:
        if not os.path.exists(p):
            os.makedirs(p)

    return data_path, result_path, index_path, search_path, eval_path, npz_path, log_path

def load_files(file_name, clean):
    data_path, result_path, index_path, search_path, eval_path, npz_path, log_path = load_paths(file_name)
    paths = (result_path, index_path, search_path, npz_path, eval_path, log_path)

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

def get_summary(file_name, options):
    data_path, result_path, index_path, search_path, eval_path, npz_path, log_path = load_paths(file_name)
    if not os.path.exists(result_path + 'num.npz'):
        seqs, vals, num_seqs, opts, Op = load_files(file_name, clean=False)
        seq_lens = [len(seqs[i]) for i in range(len(seqs))]
        np.savez(result_path+'num.npz', N=len(seqs), num_seqs = num_seqs, options = options, Op = Op, seq_lens = seq_lens)
    # check if the main parameters are consistent with the previous run
    # if target is clean, it doesn't matter, 
    # if it's continue, options will be replaced by the previous ones
    elif options.target!='clean' and options.target!='continue':
        summary = np.load(result_path + 'num.npz')
        opts = summary['options']
        opts = opts[()]
        if options.file_name != opts.file_name:
            raise(Exception('file_name incosistent: ', options.file_name, opts.file_name ))
        if options.proj_dim!=opts.proj_dim:
            raise(Exception('proj_dim is incosistent: ', options.proj_dim, opts.proj_dim));
        if options.big_k!= opts.big_k:
            raise(Exception('big_k is incosistent : ', options.big_k, opts.big_k))
        if options.small_k!= opts.small_k:
            raise(Exception('small_k is incosistent : ', options.small_k, opts.small_k))
        if options.metric!= opts.metric:
            raise(Exception('metric is incosistent : ', options.metric, opts.metric))
    

def get_job_path(job, job_paths, ext):
    if isinstance(job_paths, str): # check if it's only one path
        path = job_paths
    else:
        path = job_paths[job[0]]
    S = path + job[0]
    for i in range(1,len(job)):
        S = S + '_' + str(job[i])
    return S + ext


# generate job-specific arguments
def gen_job_args(job):
    S = ' --target ' + job[0] 
    argnames = ['', ' --seq-id ', ' --search-seq-id '] 
    for i in range(1,len(job)):
        S = S + argnames[i] + str(job[i])
    return S 



# generate the full bsub command based on job and resources 
# resources default to options defaults if not provided
def get_bsub_cmd(job, options, mem = -1, t = -1):
    _, _, _, _, _, _, log_path = load_paths(file_name)
    # if not override, use options defaults
    if mem==-1:
        mem = options.memory
    if t==-1:
        t = options.time

    # create bsub command 
    command = 'bsub -o  ' + log_path                                # save lsf.o output in the log
    command = command + ' -R "rusage[mem=' + str(mem) + '000]" '    # allocate memory 
    command = command + ' -W ' + str(t) +  ':00 '                   # allocate time 
    command = command + ' python main.py '                # name of the process

    # add options, but remove args irrelevant to worker processes
    for k,v in options.__dict__.items():
        if k not in ['target', 'seq_id', 'search_seq_id', 'target', 'forward_target', 'memory', 'time', 'directory']:
            command = command + ' --' + k.replace('_','-') + ' ' + str(v)       
    # get argname from varname 
    command = command + ' --time ' + str(t) + ' --memory ' + str(mem) + ' --directory ' + RESULT_DIR

    # add job args to bsub and option args
    command = command + gen_job_args(job)      

    return command


# recursively compute jobs that need to run based on dependency map
def crawl_dep_recurse(Map, Set, curr):
    if curr not in Set:
        Set.add(curr)
        if curr in Map.keys():
            for job in Map[curr]:
                crawl_dep_recurse(Map, Set, job)

def crawl_dep(Map, curr):
    Set = set()
    crawl_dep_recurse(Map, Set,curr)
    return Set

def add_job(Map, job, deps):
    Map[job] = deps

def run_jobs(jobs_flat, job_dep, options, job_paths):
    _, _, _, _, _, _, log_path = load_paths(file_name)
    started = dict()
    for j in jobs_flat:
        started[j] = False
    printed = dict()
    for j in jobs_flat:
        printed[j] = False
    fcommands = open(log_path + 'commands.sh', 'a+')
    flog = open(log_path + 'log.txt', 'a+')
    counter = -1 
    finished = [False] * len(jobs_flat)
    while not sum(finished)==len(finished):
        counter = (counter + 1) % 5
        # find and hanlde killed jobs (rerun with more resource)
        if counter == 0:
            killed_jobs = find_killed_jobs(log_path, flog)
            for kj in killed_jobs:
                if kj.job in jobs_flat:
                    bsub_cmd = newcmd_for_killed(fname = kj.fname, job = kj.job, 
                            memory = kj.memory, t = kj.time, reason = kj.reason, 
                            log_path = log_path, options = options)
                    bsub_out = subprocess.check_output(bsub_cmd, shell=True)
                    if options.verbose>=1:
                        print('RUNNING after killed : ', bsub_cmd)
                        sys.stdout.flush()
                        
                    print('RUNNING after killed : ', bsub_cmd, file = flog)
                    started[kj.job] = True
                    print(bsub_out, file = flog)
                    print(bsub_cmd, file = fcommands)
        for jobi, job in enumerate(jobs_flat):
            done_path = get_job_path(job,job_paths, '.done')
            outp = get_job_path(job,job_paths, '.txt')
            finished[jobi] = os.path.exists(done_path)
            # if the dependencies of "job" are satisfied dispatch it
            if not os.path.exists(done_path) and not started[job]:
                deps_satisfied = [os.path.exists(get_job_path(j,job_paths, '.done')) for j in job_dep[job] ] 
                if all(deps_satisfied):
                    bsub_cmd = get_bsub_cmd(job, options)  
                    bsub_out = subprocess.check_output(bsub_cmd, shell=True)
                    if options.verbose>=1:
                        print('RUNNING: ', bsub_cmd)
                        sys.stdout.flush()
                    print('RUNNING: ', bsub_cmd, file = flog)
                    started[job] = True
                    print(bsub_out, file = flog)
                    print(bsub_cmd, file = fcommands)
            # print jone jobs based on verbosity
            # vero verbose=2 print all outputs, for verbose = 1 only print the final result
            if options.verbose==0:
                os.system('clear')
                print('Finished jobs: ', sum(finished), '/', len(finished))
                Len = 80
                Char = '\u25A0'
                prog = Len* sum(finished) * 1.0/len(finished)
                prog = int(prog)
                print('Progress bar:')
                print('|'+Char*prog + ' '*(Len-prog) + '|')
            if (options.verbose == 2 or (options.verbose==1 and job[0]=='merge') ) and (not printed[job] and os.path.exists(outp) and os.path.exists(done_path) ):
                    f = open(outp,'r')
                    print(f.read())
                    sys.stdout.flush()
                    printed[job] = True
                    f.close()

        time.sleep(1)
    fcommands.close()
    flog.close()
    print('result_path: ', result_path, '\n')

    

# construct job tupple from input and return job
# if the job is inconsistent with job format return error=True
def get_job_from_args(target, i, j):
    if i>=0 and j>=0:
        job =  (target, i, j)
    elif i>=0:
        job =  (target, i)
    else:
        job =  (target,)
    error = False
    if target=='':
        error = True
    if target in ['search', 'eval'] and j==-1:
        error = True
    if target in ['build', 'search', 'eval'] and i==-1:
        error = True
    
    return job, error

def get_job_from_target(target, i, j):
    if target in ['clean', 'merge', 'final', 'all']:
        job = (target,)
    elif target in ['search', 'eval']:
        job = (target, i, j)
    elif target in ['build']:
        job = (target, i)
    else:
        raise(Exception('get_job_from_target used with the wrong target ', target))
    return job


# create new bsub command for the killed job 
def newcmd_for_killed(fname, job, memory, t, reason, log_path, options):
    if reason=='memory':
        cmd = get_bsub_cmd(job = job, options = options, mem = memory*2, t = t ) 
    elif reason=='time':
        if t<=4:
            t2 = 24
        elif t<=24:
            t2 = 120
        else:
            raise(Exception('the killed jobs has taken more than 120 hours, can\'t create the new job'))
        cmd = get_bsub_cmd(job = job, options = options, mem = memory, t = t2 ) 
    else: 
        raise(Exception('uknown reason for the killed job, can\'t create the new job'))

    return cmd
    

# get killed job from lsf.o* files in the log path and moves them to log_path/killed
# throws an exception if an error happens in the job recovery
def find_killed_jobs(log_path, flog):
    if not os.path.exists(log_path + 'killed'):
        os.makedirs(log_path + 'killed')
    # search for files with patern lsf* or *.out
    files = glob.glob(log_path + 'lsf*')
    files.extend(glob.glob(log_path + '*.out'))
    killed_jobs = []
    error = False
    e_message = ''
    for fpath in files:
        fname = fpath.split('/')[-1]
        with open(fpath, 'r') as f:
            content = f.read()
            if 'job killed' in content or 'Exited with exit code' in content:
                L = content.split()
                arg_vals = ['', 10, 4, -1, -1]
                arg_names = ['--'+s for s in ['target', 'memory', 'time', 'seq-id', 'search-seq-id'] ]
                __arg_names = ['-'+s for s in ['t', 'M', 'T', 'i', 'j'] ]
                for argi,arg in enumerate(arg_names):
                    if arg in L:
                        arg_vals[argi] = L[L.index(arg)+1]
                    elif __arg_names[argi] in L:
                        arg_vals[argi] = L[L.index(__arg_names[argi]) + 1]
                    if argi>0 and isinstance(arg_vals[argi], str):
                        s = re.sub(r'\W+', '', arg_vals[argi])
                        arg_vals[argi] = int(s)
                target, memory, time, i, j = arg_vals
                reason = 'unkown'
                if 'TERM_MEMLIMIT' in content:
                    reason = 'memory'
                elif 'TERM_RUNLIMIT' in content:
                    reason = 'time' 
                elif 'TERM_OWNER' in content:
                    reason = 'user'
                elif 'Exited with exit code' in content:
                    reason = 'error'
                job, e =  get_job_from_args(target, i, j)
                if e==True or reason in ['unkonwn' ,'error']:
                    error = True
                e_message = ''
                if e==True:
                    e_message = ("killed job recovered from LSF file incorrectly: "+ str(job))
                if reason == 'unkown':
                    e_message = e_message + '\n' + ("job killed for unkown reason "+ str(job) )
                if reason == 'error':
                    err_index=  content.index('The output (if any) follows:')
                    excerpt = content[err_index:]
                    e_message = e_message + '\n' + excerpt
                Dict = {'fname' : fname, 'job' : job, 'memory':memory, 'time' : time, 'reason': reason}
                kj = AttrDict(Dict)
                shutil.move( fpath, log_path + 'killed/' + fname)
                print('found killed job : ', job, ' mem = ', memory, ' time = ', time, ', and moved lsf file: ', fname, ' to killed/ ' )
                print('found killed job : ', job, ' mem = ', memory, ' time = ', time, ', and moved lsf file: ', fname, ' to killed/ ' , file = flog)
                # add job if it's not killed by the user
                if reason != 'user':
                    killed_jobs.append(kj)
    # raise an exception of an error occured in processing killed jobs
    if error == True:
        raise(Exception(e_message)) 

    return killed_jobs





if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-f','--file-name',dest='file_name', default = '', type='string', help = 'seq file to load')
    parser.add_option('-i','--seq-id',dest='seq_id', default = 0, type='int', help = 'seq id to build index for')
    parser.add_option('-j','--search-seq-id',dest='search_seq_id', default = 0, type='int', help = 'search sequence id')
    parser.add_option('-k','--small-k',dest='small_k', default = 3, type='int', help = 'small k to do kmer counting')
    parser.add_option('-K','--big-k',dest='big_k', default = 100, type='int', help = 'big K to slide over the strings')
    parser.add_option('-D','--proj-dim',dest='proj_dim', default = 120, type='int', help = 'number of dimensions to project onto')
    parser.add_option('-n','--num-trees',dest='num_trees', default = 2, type='int', help = 'number of trees to be built')
    parser.add_option('-s','--step-build',dest='step_build', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-S','--step-search',dest='step_search', default = 5, type='int', help = 'step for indices of KNN built')
    parser.add_option('-m','--metric',dest='metric', default = 'euclidean', type='string', help = 'metric to be used')
    parser.add_option('-N','--num-neighbors',dest='num_neighbors', default = 50, type='int', help = 'number of nearest neighbors to search')
    parser.add_option('-t','--target',dest='target', default = 'build', type='string', help = 'metric to be used')
    parser.add_option('-F','--forward-target',dest='forward_target', default = 'merge', type='string', help = 'metric to be used')
    parser.add_option('-M','--memory',dest='memory', default = 30, type='int', help = 'memory to allocate for each process')
    parser.add_option('-T','--time',dest='time', default = 4, type='int', help = 'designated hours for each process ')
    parser.add_option('-v','--verbose',dest='verbose', default = 2, type='int', help = ' level of verbosity (0,1, & 2) ')
    parser.add_option('-d','--directory',dest='directory', default = 'exp', type='str', help = ' directory for results (exp: save in experimental, auto: create new dir) ')
    (options, args) = parser.parse_args()

    # set the result directory based on options.directory
    if options.directory=='auto':
        RESULT_DIR = options.file_name + '_' + str(datetime.now()).replace('-','_').replace(':','_').replace(' ','_').replace('.','_')
    elif options.directory=='exp':
        RESULT_DIR = 'experimental/' + options.file_name
    else:
        RESULT_DIR = options.directory

    data_path, result_path, index_path, search_path, eval_path, npz_path, log_path = load_paths(options.file_name)
    # if target is to continue do the all target, but load the options from the file
    target = options.target
    if options.target=='continue':
        if options.directory=='auto':
            raise(Exception('target continue cannot work with automatic directory'))
        if not os.path.exists(result_path + 'num.npz'):
            raise(Exception('num.npz not found, target continue requires it to run'))
        if options.file_name == '':
            raise(Exception('target continue requires specifying --file-name '))
        summary = np.load(result_path+'num.npz')
        saved_options = summary['options']
        options = saved_options[()]
        target = 'all'

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

    job_paths = {'clean':result_path, 
            'build': index_path,
            'search' : search_path,
            'eval': eval_path,
            'merge': result_path ,
            'final': result_path , 
            'all': result_path , 
            'npz' : npz_path }

    get_summary(file_name, options)

    job_done = get_job_from_target(target, sid, sid2)
    fout_path = get_job_path(job_done, job_paths, '.txt')
    fout = open(fout_path, 'w+')


    if target=='all':
        summary = np.load(result_path+'num.npz')
        N = summary['N']

        jobs_dep = dict()

        # add clean job if directory is not set automatically
        if options.directory!='auto':
            add_job(jobs_dep,('clean',),[])
        build_jobs = []
        for i in range(N):
            job = ('build',i) 
            # add clean dependency if directory is not set automatically
            if options.directory=='auto':
                add_job(jobs_dep,job,[])
            else:
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
        run_jobs(jobs_flat, jobs_dep, options, job_paths) 

    elif target=='final':
        pass

        
    elif target=='merge':
        quads = set()
        summary = np.load(result_path+'num.npz')
        N = summary['N']
        num_seqs = summary['num_seqs']
        Op = summary['Op']
        Op = Op[()]
        total_count = 0
        total_correct = 0
        for i in range(N):
            for j in range(N):
                if i!=j:
                    res = np.load(eval_path +str(i)+'_'+str(j)+'.npz') 
                    q = res['quadrupples']
                    quads = quads.union(q[()])
                    total_count = total_count + res['total_count']
                    total_correct = total_correct + res['total_correct']
        Total = 0
        for ns in num_seqs:
            Total = Total + ns*(ns-1)/2 * Op.num_genes
        if 'seq_lens' in summary.files:
            seq_lens = summary['seq_lens']
        else:
            seqs, _, _, _, _ = load_files(file_name, clean=False)
            seq_lens = [len(seqs[i]) for i in range(len(seqs))]

        print('\nMERGE RESULT ' + '-'*50, file=fout)
        print('number of seqs: ', N , file=fout)
        print('seq options : ', Op, file = fout)
        print('analysis options : ', options, file=fout)
        print('mean seq len: ', np.mean(seq_lens), file = fout)
        print('false positive : ', (total_count - total_correct ) *1.0/total_count, file = fout)
        print('final recall : ', len(quads)*1.0/Total, file=fout)
        print('RESULT_DIR: ', RESULT_DIR, file=fout)
        print('-'*63 + '\n', file = fout)


    elif target=='clean':
        load_files(file_name, clean=True)
        get_summary(file_name, options)


    elif target=='build':
    
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

    elif target=='search':

        data1 = np.load(npz_path+'data'+str(sid)+'.npz')
        knn_index = AnnoyIndex(proj_dim, metric)
        knn_index.load(index_path + str(sid) + '.ann', prefault=True)
        build_indices = data1['build_indices']

        data2 = np.load(npz_path+'data'+str(sid2)+'.npz')
        convolved = data2['convolved']
        search_indices = np.arange(0,convolved.shape[0],step_search)
        NN, NN_dist = utility.knn_search_value(knn_index, convolved, search_indices, build_indices,num_neighbors)
        np.savez(search_path + str(sid)+'_'+str(sid2) + '.npz', NN=NN, NN_dist=NN_dist, search_indices = search_indices)

    elif target=='eval':

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

        print('\nEVAL RESULT ' + '-'*50, file=fout)
        print('eval job: -i ', sid, ', -j ', sid2, file=fout)
        print('false positive = ', (total_count - total_correct)/total_count, file=fout)
        print('recall = ', len(quadrupples)*1.0/Op.num_genes, file=fout)
        print('RESULT_DIR: ', RESULT_DIR, file=fout)
        print( '-'*62 + '\n', file=fout)
        
    # print directory of results in case it was missed 
    # print done for all targets 
    # except for all target that is the dispatcher 
    fout.close()
    if os.path.exists(fout_path) and os.path.getsize(fout_path)==0:
        os.remove(fout_path)
    Job_path = get_job_path(job_done, job_paths, '.done')
    print('dummy', file=open(Job_path, 'w+'))
