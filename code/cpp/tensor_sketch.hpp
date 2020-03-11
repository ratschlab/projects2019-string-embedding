
#include "string_options.hpp"
#include "string_tools.hpp"
#include "st_utils.hpp"


template<typename C1, typename C2> 
struct tensor_sketch : public string_tools_t {

    Vec3D<C1> rand_phase;
    Vec2D<double> distribution;


    void init_rand(int sig_len, int dim, int t_len, int num_bins) {
        std::cauchy_distribution<double> cauchy(0,1);
        std::uniform_int_distribution<C1> unif(0,num_bins-1);
        int tsize = pow(sig_len,t_len);
        distribution = Vec2D<C2>(dim, Vec<C2>(num_bins));
        for (int d=0; d<dim; d++) {
            for (int di=0; di<num_bins; di++) {
                distribution[d][di] = cauchy(gen); 
            }
        }
        rand_phase = Vec3D<C1>(dim, Vec2D<C1>(t_len, Vec<C1>(sig_len)));
        for (size_t d=0; d<dim; d++) {
            for (size_t ti=0; ti<t_len; ti++) {
                for (size_t c = 0 ; c<sig_len; c++ ) {
                    rand_phase[d][ti][c] = unif(gen);
                }
            }
        }
    }

    template <typename T1, typename T2>
    void sketch(const Vec<T1> &seq, Vec<T2> &H, int sig_len, int dim, int t_len, int num_bins, bool normalize) {
        assert(rand_phase.size()==dim);
        assert(distribution.size()==dim);
        assert(rand_phase[0].size()==t_len);
        assert(distribution[0].size()==num_bins);
        Vec3D<C1> mem(dim,Vec2D<C1>(t_len, Vec<C1>(num_bins,0)));
        for (size_t d=0; d<dim; d++ ) { 
            for (size_t i=0; i<seq.size(); i++)  {
                for (size_t t=0; t<t_len; t++) {
                        auto  ph = rand_phase[d][t][seq[i]];
                    if (t==0) {
                        mem[d][t][ph]++;
                    } else {
                        for (size_t p=0; p<num_bins; p++) {
                            assert(t<mem[d].size() and t>=1);
                            mem[d][t][(p+ph) % num_bins] += mem[d][t-1][p];
                        }
                    }
                }
            } 
        }
        H = Vec<T2>(dim, 0);
        for (int m=0; m<dim; m++) {
            for (int d=0; d<num_bins; d++) {
                H[m] += distribution[m][d]*mem[m][t_len-1][d];
            }
            H[m] = H[m]/dim;
//            H[m] = mem[m][t_len-1][0];
        }
        if (normalize) {
            for (int m=0; m<dim; m++) {
                H[m] = H[m]/L1(mem[m][t_len-1]);
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
