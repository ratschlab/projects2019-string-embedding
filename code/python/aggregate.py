import numpy as np
# import utility
from python import utility

PROJ_DIR = utility.proj_dir()


def get_quads(n, p):
    quads = set()
    for i in range(n):
        for j in range(n):
            if i!=j:
                ev = np.load(PROJ_DIR + '/results/' + p + '/eval/' + str(i) + '_' + str(j) + '.npz')
                qs = ev['quadrupples']
                quads = quads.union(qs[()])
    return quads


def get_acc(paths):
    quads = set()
    for p in paths:
        q = get_quads(5,p)
        quads= quads.union(q)
    return len(quads)/60000


if __name__ == '__main__':
    paths = ['seqs12_2019_06_16_05_51_11_641954', 'seqs12_2019_06_16_05_50_46_528370', 'seqs12_2019_06_16_05_50_46_109711','seqs12_2019_06_16_05_50_12_881948', 'seqs12_2019_06_15_14_29_46_040560', 'seqs12_2019_06_15_14_29_46_040549'] 
    print(' accuracy for seq12, dim=10 = ', get_acc(paths))
