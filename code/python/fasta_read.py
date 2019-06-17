import numpy as np
import glob
from attrdict import AttrDict
import utility
import sys

def load_seqs_vals(file_name, i):
    # read maf file 
    proj_dir = utility.proj_dir() + '/data/fasta/'
    datadir = proj_dir + file_name
    fasta_files = glob.glob(datadir+'/seqs*.fa')
    N = len(fasta_files)
    genes = [[] for i in range(N)]
    maf_file = open(datadir + '/MSA.maf','r')
    line = maf_file.readline()
    while '===' not in line:
        if line[0]=='s':
            line = maf_file.readline()
            continue
        seq, g, s, l = [int(i) for i in line.split()]
        genes[seq].append([g, s, l])
        line = maf_file.readline()

    f = open(fasta_files[i])
    f.readline()
    seq = f.readline()[:-1]
    L = len(seq)
    val = np.zeros(L, dtype=np.uint32)
    for gene in genes[i]:
        gc, s, l = gene
        val[s:s+l] = gc

    return seq, val


def load_options(file_name):
    proj_dir = utility.proj_dir() + '/data/fasta/'
    datadir = proj_dir + file_name
    opts_files = open(datadir+'/options.txt')
    lines = opts_files.readlines()
    for line in lines:
        if 'num_seqs' in line:
            num_seqs = [int(i) for i in line.split()[1:]]
        elif 'gene_lens' in line:
            gene_lens = [int(i) for i in line.split()[1:]]
        elif 'seq_lens' in line:
            seq_lens = [int(i) for i in line.split()[1:]]
        elif 'seq_options' in line:
            opts = line.replace(',', ' ').replace('}','').replace('{','').split()[2:]
            Dict = dict()
            for i in range(int(len(opts)/2)):
                key, val = opts[2*i:2*i +2]
                if key in ['mutation_rate', 'geometric_p', 'reverse_p']:
                    Dict[key] = float(val)
                elif key == 'file_name':
                    Dict['save_directory'] = val
                else:
                    Dict[key] = int(val)
    Dict['num_seqs'] = num_seqs
    Dict['gene_lens'] = gene_lens
    Dict['seq_lens'] = seq_lens
    Dict['N'] = len(seq_lens)
    options = AttrDict(Dict)    
    return Dict


if __name__ == '__main__':
    print(sys.argv)

    file_name = sys.argv[1]
    print('file_name = ', file_name)
    proj_dir = utility.proj_dir() + '/data/fasta/'
    datadir = proj_dir + file_name

    # read maf file 
    fasta_files = glob.glob(datadir+'/seqs*.fa')
    N = len(fasta_files)
    genes = [[] for i in range(N)]
    maf_file = open(datadir + '/MSA.maf','r')
    line = maf_file.readline()
    while '===' not in line:
        if line[0]=='s':
            line = maf_file.readline()
            continue
        seq, g, s, l = [int(i) for i in line.split()]
        genes[seq].append([g, s, l])
        line = maf_file.readline()

    # read fasta files and fill the values 
    seqs = []
    vals = []

    for i in range(N):
        f = open(fasta_files[i])
        f.readline()
        seq = f.readline()[:-1]
        L = len(seq)
        val = np.zeros(L, dtype=np.uint32)
        for gene in genes[i]:
            gc, s, l = gene
            val[s:s+l] = gc
        seqs.append(seq)
        vals.append(val)



    # read options file 

    opts_files = open(datadir+'/options.txt')
    lines = opts_files.readlines()
    for line in lines:
        if 'num_seqs' in line:
            num_seqs = [int(i) for i in line.split()[1:]]
        elif 'gene_lens' in line:
            gene_lens = [int(i) for i in line.split()[1:]]
        elif 'seq_options' in line:
            opts = line.replace(',', ' ').replace('}','').replace('{','').split()[2:]
            Dict = dict()
            for i in range(int(len(opts)/2)):
                key, val = opts[2*i:2*i +2]
                if key in ['mutation_rate', 'geometric_p', 'reverse_p']:
                    Dict[key] = float(val)
                elif key == 'file_name':
                    Dict['save_directory'] = val
                else:
                    Dict[key] = int(val)
    Dict['num_seqs'] = num_seqs
    Dict['gene_lens'] = gene_lens
    options = AttrDict(Dict)    


    # save the values
    np.savez(utility.proj_dir() + '/data/'+file_name+'.npz', seqs = seqs, vals = vals, num_seqs = num_seqs, gene_lens = gene_lens, options = Dict )
