
#include <numeric>

#include "string_options.hpp"
#include "string_tools.hpp"


template<typename C1, typename C2> 
struct tensor_embed : public string_tools_t {

    Vec2D<C2> proj_ideal;

    void init_ideal(int sig_len, int t_len, int dim) {
        int tsize = pow(sig_len,t_len);
        proj_ideal = Vec2D<C2>(dim, Vec<C2>(tsize));
        std::cauchy_distribution<C2> cauchy(0,1);
        for (int d=0; d<dim; d++) {
            for (int t=0; t<tsize; t++) {
                proj_ideal[d][t] = cauchy(gen);
            }
        }
    }

    template <typename T1, typename T2>
    void embed_naive(const Vec<T1> &seq, Vec<T2>  &T, int sig_len, int t_len) {
        int tsize = pow(sig_len,t_len);
        T = Vec<T2>(tsize,0);
        Vec<size_t> ind(t_len); 
        std::iota(ind.begin(), ind.end(), 0);
        do {
            auto ti = read_tuple(seq,ind,sig_len);
            assert(ti<T.size());
            T[ti]++;
        } while (inc_ind_sorted(ind, seq.size()));
    }

    template <typename T1, typename T2> 
    void sketch_ideal(const Vec<T1> &T, Vec<T2> &h, int num_bins, bool normalize) {
        int dim = proj_ideal.size();
        int tsize = proj_ideal[0].size();
        h = Vec<T2>(dim, 0);
        for (size_t d=0; d<dim; d++) {
            for (size_t i=0; i<T.size() ; i++ ) {
                h[d] += proj_ideal[d][i]*T[i] ;
            }
            h[d] = h[d]/dim;
        }
        if (normalize) {
            for (size_t d=0; d<dim; d++ ) { 
                h[d] = h[d] / L1(T);
                h[d] = h[d] * num_bins;
            }
        }
    }

    string_opts ops;
    explicit tensor_embed(const string_opts &ops);


    template<typename T1, typename T2>
    void sketch_ideal(const Vec <T1> &seq, Vec <T2> &H) {
        sketch_ideal(seq, H, ops.num_bins, ops.normalize);
    }

    void init_ideal() {
        init_ideal(ops.sig_len, ops.t_len, ops.dim);
    }

    template<typename T1, typename T2>
    void embed_naive(const Vec <T1> &seq, Vec <T2> &T) {
        embed_naive(seq, T, ops.sig_len, ops.t_len);
    }
};

template<typename C1, typename C2>
tensor_embed<C1, C2>::tensor_embed(const string_opts &ops) : ops(ops) {}
