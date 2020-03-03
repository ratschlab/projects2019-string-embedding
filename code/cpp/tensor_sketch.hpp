
#include "st_utils.hpp"
#include "string_options.hpp"
#include "string_tools.hpp"


template<typename C1, typename C2> 
struct tensor_sketch : public string_tools_t {

    Vec2D<C1> rand_phase;
    Vec2D<C2> distribution;


    void init_rand(int sig_len, int dim, int t_len, int num_bins) {
        std::cauchy_distribution<double> cauchy(0,1);
        std::uniform_int_distribution<int> unif(0,num_bins-1);
        int tsize = sutils.pow(sig_len,t_len);
        rand_phase = Vec2D<C1>(t_len, Vec<C1>(sig_len));
        distribution = Vec2D<C2>(dim, Vec<C2>(num_bins));
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

    template <typename T1, typename T2>
    void sketch(const Vec<T1> &seq, Vec<T2> &H, int sig_len, int dim, int t_len, int num_bins, bool normalize) {
        assert(rand_phase.size()==t_len);
        assert(distribution.size()==dim);
        assert(rand_phase[0].size()==sig_len);
        assert(distribution[0].size()==num_bins);
        Vec2D<C1> mem(t_len, Vec<C1>(num_bins,0));
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
        H = Vec<T2>(dim, 0);
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
                H[m] = H[m]/sutils.L1(mem[t_len-1]);
                H[m] = H[m]/num_bins;
            }
        }
    }


    string_opts ops;
    tensor_sketch(const string_opts ops): ops(ops)  {}


    template <typename T1, typename T2>
    void sketch(const Vec<T1> &seq, Vec<T2> &H ) {
        sketch(seq, H, ops.sig_len, ops.dim, ops.t_len, ops.num_bins, ops.normalize);
    }

    void init_rand() {
        init_rand(ops.sig_len, ops.dim, ops.t_len, ops.num_bins);
    }
};
