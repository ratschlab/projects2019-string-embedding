import random
import sys
from optparse import OptionParser
import numpy as np

# colors to be used in the terminal
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
colors = [bcolors.HEADER, bcolors.OKBLUE, bcolors.OKGREEN, bcolors.WARNING, bcolors.FAIL]


# genrate random sequence 
def random_seq(ref,L):
    S = [random.choice(ref) for i in range(L)]
    return ''.join(S)



# mutate sequence gien the mutatioin rate p
def mutate_seq2(seq, val , ref = 'acgt', p = 0.01, ps = [1/3, 1/3, 1/3], reverse_p = .1, geometric_p = 1):
    seq2 = ''
    val2 = []
    i = 0
    while (i < len(seq)):
        r = random.random()
        if r<p:
            j = np.nonzero(np.cumsum(ps)>random.random())[0][0]
            #j = random.choice([0, 1, 2])
            
            rounds = np.random.geometric(geometric_p)
            for r in range(rounds):
                if i==len(seq):
                    break
                # change
                if j==0:
                    seq2 += random.choice(ref)
                    if seq[-1]==seq[i]:
                        val2.append(val[i])
                    else:
                        val2.append(0)
                    i = i + 1
                # insert
                elif j==1:
                    seq2 += random.choice(ref)
                    val2.append(0)
                    # delete
                elif j==2:
                    i = i + 1
                    pass
        else:
            seq2 += seq[i]
            val2.append(val[i])
            i = i + 1
    if random.random()<reverse_p:
        seq2 = seq2[::-1]
        val2 = val2[::-1]
    return seq2, val2



# generate high variation sequence given the parameters
def high_var_seq2(num_genes, padding, gene_len, num_seq, mutation_rate, repeat, 
                 gene_len2 = -1, 
                 num_seq2 = -1,
                 geometric_p = 1,
                 reverse_p = 0, 
                 print_seqs = False, colored = False, print_vals = False):
    def gene_code(v):
        if v==0:
            return ' '
        else:
            return str(v)
    if gene_len2<0:
        gene_len2 = gene_len
    if num_seq2<0:
        num_seq2 = num_seq
    
    gene_lens = np.random.randint(gene_len,gene_len2+1,size=num_genes)
    num_seqs = np.random.randint(num_seq,num_seq2+1,size=repeat)

    seqs = []
    vals = []

    for ri in range(repeat):
        gene = []
        for k in range(num_genes):
            gene .append(random_seq('acgt',gene_lens[k]))
        for i in range(num_seqs[ri]):
            seq = ''
            val = []
            for k in range(num_genes):
                seq += random_seq('acgt',padding)
                val.extend([0]*padding)
                seq += gene[k]
                val.extend([(ri*num_genes)+k+1]*gene_lens[k])
                seq += random_seq('acgt',padding)
                val.extend([0]*padding)
            if i>0:
                seq, val = mutate_seq2(seq,val, p=mutation_rate, reverse_p = reverse_p, geometric_p = geometric_p)
            seqs.append(seq)
            vals.append(val)

    if print_seqs == True:
        for i in range(len(seqs)):
            for c, v in zip(seqs[i], vals[i]):
                if v>0 and colored:
                    sys.stdout.write(bcolors.UNDERLINE + colors[v % 5] + c + bcolors.ENDC)

                else:
                    sys.stdout.write(c)
            print("")
            if print_vals:
                print(''.join([gene_code(v) for v in vals[i]]))
                print("")

    return seqs, vals, gene_lens, num_seqs 




# sample command line runs : python high_var_seqs.py  --mutation-rate .1  --geometric-p .4 --num-seq 2 --num-seq2 5  --gene-len 100 --gene-len2 400 -G 50 --repeat 3 --padding 1000 --save-directory ../data/seqs7

if __name__ == '__main__':
    # reading the arguments
    #print('generating the sequences')
    parser = OptionParser()
    parser.add_option('-m','--mutation-rate',dest='mutation_rate', default = .2, type='float', help = 'mutation rate per basepair, range: (0,1)')
    parser.add_option('-g','--gene-len',dest='gene_len', default = 10, type='int', help = 'length of each gene segment')
    parser.add_option('-G','--num-genes',dest='num_genes', default = 1, type='int', help = 'number of genes per genome (per sequence)')
    parser.add_option('-P','--padding',dest='padding', default = 200, type='int', help = 'number of iid random basepairs between gene segments')
    parser.add_option('-n','--num-seq',dest='num_seq', default = 2, type='int', help =  'number of genomes/sequences')
    parser.add_option('-r','--repeat',dest='repeat', default = 1, type='int', help = 'number of iid dupplicates of the whole datasets')
    parser.add_option('-Z','--geometric-p',dest='geometric_p', default = 1, type='float', help = 'l factor for geometric distribution for indels')
    parser.add_option('-V','--reverse-p',dest='reverse_p', default = 0, type='float', help = 'probability for reversing a gene')
    parser.add_option('-X','--gene-len2',dest='gene_len2', default = -1, type='int', help = 'maximum gene length, defaults to gene-len')
    parser.add_option('-Y','--num-seq2',dest='num_seq2', default = -1, type='int', help = 'maximum number of sequence in each repeat block, defaults to num-seq')
    parser.add_option('-c','--colored',dest='colored', default = False,action='store_true', help = 'color basepairs corresponding to gene segments (mod 5)')
    parser.add_option('-v','--print-values',dest='print_values', default = False,action='store_true', help = 'print values corresponding to gene number')
    parser.add_option('-p','--print-seqs',dest='print_seqs', default = False,action='store_true', help = 'print sequences values')
    parser.add_option('-d','--save-directory',dest='save_directory', default = '../data/synth_seqs', type='string', help = 'directory to save gene')

    (options, args) = parser.parse_args()

    seqs, vals, gene_lens, num_seqs = high_var_seq2(options.num_genes, options.padding, options.gene_len, 
                                                    options.num_seq, options.mutation_rate, 
                                                    repeat = options.repeat, 
                                                    colored = options.colored, 
                                                    print_seqs = options.print_seqs, 
                                                    print_vals=options.print_values, 
                                                    gene_len2 = options.gene_len2, 
                                                    num_seq2 = options.num_seq2,
                                                    geometric_p = options.geometric_p, 
                                                    reverse_p = options.reverse_p )
    
    np.savez(options.save_directory, seqs=seqs, vals=vals, options=options, num_seqs = num_seqs, gene_lens = gene_lens)
    #high_var_seq(options.num_genes, options.padding, options.gene_len, options.num_seq, options.mutation_rate, repeat = options.repeat, colored = options.colored, print_seqs = options.print_seqs, print_vals=options.print_values)



