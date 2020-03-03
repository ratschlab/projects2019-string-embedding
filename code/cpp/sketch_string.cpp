#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include "command_group.h"
#include "st_utils.hpp"

#include "string_tools.cpp"


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


struct test_unit { 
        typedef tensor_embed<int,int> TE;
        typedef string_tools_t<int,int> ST;
        typedef TE::seq_t seq_t;
        typedef TE::ivect_t ivect_t;
        typedef TE::dvect_t dvect_t;
    void test_pairs(std::string alpha, int num, string_opts &opts) {
        st_utils sut;
        std::ofstream  fout(opts.res_path), fconf(opts.conf_path);
        std::string sys_cmd = "mkdir -p " + opts.dir;
        std::system(sys_cmd.data());
        fconf << opts.get_config(1);
        std::cout << opts.get_config(0) ;

        std::vector<seq_t> S1, S2;
        ST ss;
        TE st(opts);
        opts.sig_len = sut.pow(opts.sig_len ,opts.k_len);
        TE tt(opts);
        
        tt.init_rand();
        tt.init_ideal();

        st.gen_pairs(S1,S2,num);

        auto header = "ed, \t td, \t hd, \t HD, \n";
        std::cout << header;
        fout << header;
        for (int i=0; i<num; i ++ ){
            seq_t r1 = S1[i], r2 = S2[i]; 
            ivect_t s1, s2, T1, T2;
            dvect_t h1, h2, H1, H2;

            st.seq2kmer(r1,s1);
            st.seq2kmer(r2,s2);

            tt.embed_naive(s1,T1);
            tt.embed_naive(s2,T2);
            tt.sketch_ideal(T1,h1);
            tt.sketch_ideal(T2,h2);
            tt.sketch(s1,H1);
            tt.sketch(s2,H2);

            auto tdiff = sut.L1_diff<int>(T1,T2);
            auto hdiff = sut.median<double>(h1,h2);
            auto Hdiff = sut.median<double>(H1,H2);
            auto ed = ss.edit_distance(r1,r2);
            std::cout << ed << ",\t " << tdiff << ",\t " <<  
                std::setprecision(4) <<  hdiff<< ",\t " << std::setprecision(4) << Hdiff<< "\n";
            fout << ed << ",\t " << tdiff << ",\t " <<  
                std::setprecision(4) <<  hdiff<< ",\t " << std::setprecision(4) << Hdiff<< "\n";
        }
    }
};


int main(int argc, char* argv[]) {
    using std::cout;
    string_opts  opts;
    test_unit tu;
    opts.read_args(argc,argv);

    int num = 100;
    std::string alpha = "acgt"; 
    tu.test_pairs(alpha, num, opts);

    return 0;
}
