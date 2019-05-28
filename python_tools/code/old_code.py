
# mutate sequence gien the mutatioin rate p
def mutate_seq(seq, val, ref = 'acgt', p = 0.01, ps = [1/3, 1/3, 1/3]):
    seq2 = ''
    val2 = []
    for i in range(len(seq)):
        r = random.random()
        if r<p:
            j = np.nonzero(np.cumsum(ps)>random.random())[0][0]
            #j = random.choice([0, 1, 2])

            # change
            if j==0:
                seq2 += random.choice(ref.replace(seq[i],''))
                val2.append(0)
            # insert
            elif j==1:
                seq2 += random.choice(ref)+seq[i]
                val2.append(0)
                val2.append(val[i])
                # delete
            elif j==2:
                pass
        else:
            seq2 += seq[i]
            val2.append(val[i])
    return seq2, val2


def mutate_seq_align(seq, val , ref = 'acgt', p = 0.01, ps = [1/3, 1/3, 1/3], reverse_prob = .1, geometric_p = 1):
    seq2 = ''
    seq1_ = ''
    seq2_ = ''
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
                    seq1_ = seq1_ + seq[i]
                    seq2_ = seq2_ + seq2[-1]
                    if seq[-1]==seq[i]:
                        val2.append(val[i])
                    else:
                        val2.append(0)
                    i = i + 1
                # insert
                elif j==1:
                    seq2 += random.choice(ref)
                    seq1_ = seq1_ + '-'
                    seq2_ = seq2_ + seq2[-1]
                    val2.append(0)
                    # delete
                elif j==2:
                    seq1_ = seq1_ + seq[i]
                    seq2_ = seq2_ + '-'
                    i = i + 1
                    pass
        else:
            seq2 += seq[i]
            seq1_ = seq1_ + seq[i]
            seq2_ = seq2_ + seq2[-1]
            val2.append(val[i])
            i = i + 1
    if random.random()<reverse_prob:
        seq2 = seq2[::-1]
        val2 = val2[::-1]
    return seq2, val2, seq1_, seq2_




# generate high variation sequence given the parameters
def high_var_seq(num_genes, padding, gene_len, num_seq, mutation_rate, repeat, print_seqs = False, colored = False, print_vals = False):
    def gene_code(v):
        if v==0:
            return ' '
        else:
            return str(v)

    seqs = []
    vals = []

    for ri in range(repeat):
        gene = []
        for k in range(num_genes):
            gene .append(random_seq('acgt',gene_len))
        for i in range(num_seq):
            seq = ''
            val = []
            for k in range(num_genes):
                seq += random_seq('acgt',padding)
                val.extend([0]*padding)
                seq += gene[k]
                val.extend([(ri*num_genes)+k+1]*gene_len)
                seq += random_seq('acgt',padding)
                val.extend([0]*padding)
            if i>0:
                seq, val = mutate_seq(seq,val, p=mutation_rate)
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

    return seqs, vals


######### ROC curve for various thresholds for score
def ROC_curve():
    # percents = 
    accs = []
    false_pos = []
    for percent in range(0,50,5):

        th = np.percentile(scores_hit,percent) 
        pairs = []
        for ki in range(len(kmer_pos)):
            kval = kmer_vals[ki]
            pos = kmer_pos[ki]
            sid = kmer_seq_id[ki]
            this_bin = int(pos/Bin)
            if scores[sid][this_bin]>th and kval>0:
                pairs.append((sid,kval))
        pairs = set(pairs)
        acc = len(pairs)/len(seqs)/Op.num_genes
        false = sum(scores_loss>th) / ( sum(scores_loss>th) + sum(scores_hit> th) )
        accs.append(acc)
        false_pos.append(false)
        print('th ' , th, ' acc ' , acc, ' false' , false)


    plt.plot(false_pos, accs) 
    plt.ylabel('accuracy')
    plt.xlabel('false positive')

    print(len(scores_loss), len(scores_hit))
    print(sum(scores_loss>th), sum(scores_hit>th))
    print(len(kmers))
    
    
####### supposed to compute hit and loss scores
def hit_loss_plot():
    NN = np.array(NN)
    NN_dist = np.array(NN_dist)
    kmer_seq_id = np.array(kmer_seq_id)
    # total_len = len(NN)
    radius = np.percentile(NN_dist.flatten(),10)

    scores = []
    Bin = 200
    for i in range(len(seqs)):
        slen = len(seqs[i])
        scores.append([0]*int(np.ceil(slen/Bin)+1))
    # scores = np.array(scores, dtype=np.float64)

    for ii in tqdm(range(NN.shape[0])):
        ki = search_indices[ii]
        nn = NN[ii]
        dists = NN_dist[ii]

        pos = kmer_pos[ki]
        sid = kmer_seq_id[ki]
        indices = np.nonzero(kmer_seq_id[nn]!=sid)[0]
        for ni,nd in zip(nn[indices],dists[indices]):
            pos = kmer_pos[ni]
            sid = kmer_seq_id[ni]
    #         print('sid: ', sid, 'pos/Bin: ',int(pos/Bin), ' scores.shape=',scores.shape)
            scores[sid][int(pos/Bin)] = scores[sid][int(pos/Bin)] + np.exp(-nd*1.0/radius)

    scores_hit = []
    scores_loss = []
    for ki in range(len(kmer_vals)):
        kval = kmer_vals[ki]
        pos = kmer_pos[ki]
        sid = kmer_seq_id[ki]
        this_bin = int(pos/Bin)
        if kval>0:
            scores_hit.append(scores[sid][this_bin])
        else:
            scores_loss.append(scores[sid][this_bin])


    scores_hit = np.array(scores_hit)
    scores_loss = np.array(scores_loss)

    plt.figure(figsize=(12,7))
    _ = plt.hist(scores_loss, 40,density=True, alpha = .7)
    _ = plt.hist(scores_hit, 40, density=True, alpha = .7)
    plt.legend(['loss', 'hit']) 