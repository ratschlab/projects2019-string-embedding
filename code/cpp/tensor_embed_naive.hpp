
#include "st_utils.hpp"
#include "string_options.hpp"
#include "string_tools.hpp"


template<typename C1, typename C2> 
struct tensor_embed_naive : public string_tools_t<C1,C2> {

    typedef typename string_tools_t<C1,C2>::seq_t  seq_t;
    typedef typename string_tools_t<C1,C2>::ivect_t  ivect_t;
    typedef typename string_tools_t<C1,C2>::dvect_t  dvect_t;
    typedef typename string_tools_t<C1,C2>::mat2d_int  mat2d_int;
    typedef typename string_tools_t<C1,C2>::mat2d_double  mat2d_double;
    typedef typename string_tools_t<C1,C2>::indices_t  indices_t;
    typedef std::string string; 
    std::default_random_engine gen;
    st_utils sut;

    mat2d_double distribution;
    mat2d_double proj_ideal;


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
                h[i] = h[i] / sut.L1(T);
                h[i] = h[i] * num_bins;
            }
        }
    }

    string_opts ops;
    tensor_embed_naive(const string_opts ops) : ops(ops) {}


    void sketch_ideal(const ivect_t &seq, dvect_t &H ) {
        sketch_ideal(seq,H, ops.num_bins, ops.normalize);
    }

    void init_ideal() {
        init_ideal(ops.sig_len, ops.t_len, ops.dim) ;
    }
    void embed_naive(const seq_t &seq, ivect_t &T) {
        embed_naive(seq, T, ops.sig_len, ops.t_len);
    }
};
