
def get_table(filename):
    names = set()                                                                                                                                                                                       
    table = []
    with open(filename,'rb') as f:                                                                                                                                                                    
        for L in f:                                                                                                                                                                                     
            S = L.decode('utf-8')                                                                                                                                                                       
            if S[0]=='s':                                                                                                                                                                               
                row = []
                SS = S.split()                                                                                                                                                                          
                names.add(SS[1].split('.')[0])  
                row = [SS[1], SS[1].split('.')[0], SS[1].split('.')[1], int(SS[2]), int(SS[3]), SS[4]=='+', SS[6].lower()]
                table[-1][2].append(row[4])
                table[-1][0].append(row)
            if S[0]=='a':
                table.append([[], float(S.split()[1].split('=')[1]), []])

        for block in table:
            block[2] = sum(block[2])/len(block[2])


    col_names = [['full_name','name','chromosome','pos','len','dir','str'], 'score', 'mean_len']
    #del table[0]

    return table, col_names 



def table_filt(table, min_len=0, score = 0, seq_num_min=2, seq_num_max = 200):
    t = []
    for i in range(len(table)):
        if table[i][1]>=abs(score) and table[i][2]>=min_len and len(table[i][0])>=seq_num_min and len(table[i][0])<=seq_num_max:
            t.append(table[i])
    return t


def table_summary(table):
    names = set()
    full_names = set()
    for t in table:
        for r in t[0]:
            full_names.add(r[0])
            names.add(r[0].split('.')[0])
    return names,full_names



def get_chromosomes(full_names, name):
    chs= set()
    for fn in full_names:
        if name in fn:
            chs.add(fn.split('.')[1])
    return chs



def seq_vals(msa_file = 'ce10/multiz9way/chrI.maf', full_names = ['ce10.chrI', 'cb4.chrI'], root_path = '/cluster/home/ajoudaki', to_lower = True, min_len = 60):
    vals = []
    seqs = []
    table, col_names = get_table(root_path + '/genomes/' + msa_file)
    blocks = table_filt(table, min_len)
    for full_name in full_names:
        name = full_name.split('.')[0]
        chromosome = full_name.split('.')[1]
        f = open(root_path  + '/genomes/' + name + '/' + chromosome + '.fa')
        Lines = f.readlines()
        if to_lower:
            seq = ''.join([l[:-1].lower() for l in Lines[1:]])
        else:
            seq = ''.join([l[:-1] for l in Lines[1:]])
        intervals = []                         
        for blockid, block in enumerate(blocks):               
            if set(full_names).issubset([rs[0] for rs in block[0]]):  
                for r in block[0]:          
                    if r[0]==full_name:
                        intervals.append((r[3], r[3]+r[4],blockid+1))  
        intervals.sort(key = lambda tup: tup[0])
        vals.append(intervals)
        seqs.append(seq)
    return seqs, vals

def kmers_vals(seqs,vals, k = 30):
    kmers_list = []
    kmer_vals = []
    kmer_ids = []
    for si in range(len(seqs)):
        seq = seqs[si]
        intervals = vals[si]
        intervals.append((1e20,1.1e20,0)) # add an infinite interval at the end
        ii = 0
        for i in range(0,len(seq)-k+1):
            shared,comp = interval_intersect(i,i+k,intervals[ii][0],intervals[ii][1])
            V = 0
            if float(shared)/k >0.8:
                V = intervals[ii][2]
            if comp==3:
                ii += 1
            kmers_list.append(seq[i:i+k])
            kmer_vals.append(V)
            kmer_ids.append(si)

        for j in range(1,k+1):
            kmers_list.append('x'*j + seq[:k-j])
            kmer_vals.append(0)
            kmer_ids.append(si)
        for j in range(1, k):
            kmers_list.append(seq[-j:] + 'x'*(k-j))
            kmer_vals.append(0)
            kmer_ids.append(si)

    id_map = {}                                                 
    for k in kmers_list:
        id_map[k] = []
    for k,id in zip(kmers_list,kmer_ids):
        id_map[k].append(id)

    return kmers_list, kmer_vals,id_map



def interval_intersect(i1,j1,i2,j2):
    assert j1>i1 and j2>i2, "intervals not valid ... " 
    if j1<=i2:
        return 0,1
    elif i1<j2:
        return min(j1,j2)-max(i1,i2),2
    elif i1>=j2:
        return 0,3

