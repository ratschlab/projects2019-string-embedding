
#include "st_utils.hpp"
#include "command_group.h"
#include "string_tools.hpp"


template<typename C1, typename C2> 
struct tensor_embed : public string_tools_t<C1,C2> {

    typedef typename string_tools_t<C1,C2>::seq_t  seq_t;
    typedef typename string_tools_t<C1,C2>::ivect_t  ivect_t;
    typedef typename string_tools_t<C1,C2>::dvect_t  dvect_t;
    typedef typename string_tools_t<C1,C2>::mat2d_int  mat2d_int;
    typedef typename string_tools_t<C1,C2>::mat2d_double  mat2d_double;
    typedef typename string_tools_t<C1,C2>::indices_t  indices_t;
    typedef std::string string; 
    std::default_random_engine gen;
    st_utils sut;

    mat2d_int rand_phase; 
    mat2d_double distribution;
    mat2d_double proj_ideal;


    void init_rand(int sig_len, int dim, int t_len, int num_bins) {
        std::cauchy_distribution<double> cauchy(0,1);
        std::uniform_int_distribution<int> unif(0,num_bins-1);
        int tsize = sut.pow(sig_len,t_len);
        rand_phase = mat2d_int(t_len, ivect_t(sig_len));
        distribution = mat2d_double(dim, dvect_t(num_bins)); 
        for (int d=0; d<dim; d++) {
            for (int di=0; di<num_bins; di++) {
                distribution[d][di] = cauchy(gen); 
            }
        }
        for (int ti=0; ti<t_len; ti++) {
            for (int c = 0 ; c<sig_len; c++ ) {
                rand_phase[ti][c] = unif(gen);
            }
        }
    }

    void sketch(const seq_t &seq, dvect_t &H, int sig_len, int dim, int t_len, int num_bins, bool normalize) {
        assert(rand_phase.size()==t_len);
        assert(distribution.size()==dim);
        assert(rand_phase[0].size()==sig_len);
        assert(distribution[0].size()==num_bins);
        mat2d_int mem(t_len, ivect_t(num_bins,0));
        for (size_t i=0; i<seq.size(); i++)  {
            for (size_t t=0; t<t_len; t++) {
                    int ph = rand_phase[t][seq[i]];
                if (t==0) {
                    mem[t][ph]++;
                } else {
                    for (size_t p=0; p<num_bins; p++) {
                        assert(t<mem.size() and t>=1);
                        mem[t][(p+ph) % num_bins] += mem[t-1][p];
                    }
                }
            }
        } 
        H = dvect_t(dim, 0);
        for (int m=0; m<dim; m++) {
            int sum = 0;
            for (int d=0; d<num_bins; d++) {
                sum += mem[t_len-1][d];
                H[m] += distribution[m][d]*mem[t_len-1][d];
            }
            H[m] = H[m]/dim;
        }
        if (normalize) {
            for (int m=0; m<dim; m++) {
                H[m] = H[m]/sut.L1<int>(mem[t_len-1]);
                H[m] = H[m]/num_bins;
            }
        }
    }

    void init_ideal(int sig_len, int t_len, int dim) {
        int tsize = pow(sig_len,t_len);
        proj_ideal = mat2d_double(dim, dvect_t(tsize));
        std::cauchy_distribution<double> cauchy(0,1);
        for (int d=0; d<dim; d++) {
            for (int t=0; t<tsize; t++) {
                proj_ideal[d][t] = cauchy(gen);
            }
        }
    }

    void embed_naive(const seq_t &seq, ivect_t  &T, int sig_len, int t_len) {
        int tsize = pow(sig_len,t_len);
        T = ivect_t(tsize,0);
        indices_t ind(t_len); 
        for (int i=0; i<t_len; i++) { 
            ind[i]=i;
        }
        do {
            auto ti = string_tools_t<C1,C2>::read_tuple(seq,ind,sig_len);
            assert(ti<T.size());
            T[ti]++;
        } while (string_tools_t<C1,C2>::inc_ind_sorted(ind, seq.size()));
    }

    void sketch_ideal(const ivect_t &T, dvect_t &h, int num_bins, bool normalize) {
        int dim = proj_ideal.size();
        int tsize = proj_ideal[0].size();
        h.clear();
        for (int i=0; i<dim; i++) {
            double r = 0;
            for (int ti=0; ti<T.size() ; ti++ ) {
                r += proj_ideal[i][ti]*T[ti] ;
            }
            h.push_back( r );
            h[i] = h[i]/dim;
        }
        if (normalize) {
            for (int i=0; i<dim; i++ ) { 
                h[i] = h[i] / sut.L1<int>(T);
                h[i] = h[i] * num_bins;
            }
        }
    }

    // convenience functions 
    string_opts ops;
    tensor_embed(const string_opts ops) : ops(ops) {}

    void gen_pairs(std::vector<seq_t> &seqs1, std::vector<seq_t> &seqs2, int num_pairs) {
       string_tools_t<C1,C2>::gen_pairs(seqs1,seqs2,ops.len, ops.sig_len, num_pairs); 
    }

    void sketch(const seq_t &seq, dvect_t &H ) {
        sketch(seq, H, ops.sig_len, ops.dim, ops.t_len, ops.num_bins, ops.normalize);
    }

    void sketch_ideal(const ivect_t &seq, dvect_t &H ) {
        sketch_ideal(seq,H, ops.num_bins, ops.normalize);
    }

    void seq2kmer(seq_t &seq, ivect_t  &kseq) { 
        string_tools_t<C1,C2>::seq2kmer(seq, kseq, ops.k_len, ops.sig_len);
    }
    void init_rand() {
        init_rand(ops.sig_len, ops.dim, ops.t_len, ops.num_bins);
    }
    void init_ideal() {
        init_ideal(ops.sig_len, ops.t_len, ops.dim) ;
    }
    void embed_naive(const seq_t &seq, ivect_t &T) {
        embed_naive(seq, T, ops.sig_len, ops.t_len);
    }
};
