#include <algorithm>
#include <queue> 
#include <utility>
#include "st_utils.hpp"
#include "string_options.hpp"
#include "string_tools.hpp"


template<typename C1, typename C2> 
struct omp_sketch : public string_tools_t {

    Vec2D<C2> rand_permute;

    void init_permute(int sig_len, int len, int dim) {
        int uniq_siglen = (sig_len * len);
        rand_permute = Vec2D<C2>(dim, Vec<C2>(uniq_siglen));
        for (int d=0; d<dim; d++) {
            for (int t=0; t<uniq_siglen; t++) {
                rand_permute[d][t] = t;
            }
            auto rng = std::default_random_engine {};
            std::shuffle(rand_permute[d].begin(), rand_permute[d].end(), rng);
        }
    }


    void sketch(const Vec<C1> &seq, Vec2D<C2>  &T, int sig_len, int t_len) {
        size_t len = seq.size();
        size_t dim = rand_permute.size();
        size_t uniq_siglen = (sig_len * len);
        assert(rand_permute[0].size()>=uniq_siglen);
//        T = Vec2D<C2>(dim, Vec<C2>());
        T = Vec2D<C2>(dim);
        for (size_t d=0; d<dim; d++) {
            Vec<size_t> counts(sig_len,0);
            Vec<std::tuple<C2, C1, size_t>> tup_vec;
            for (size_t i=0; i<len; i++) {
                C1 Char = seq[i]; 
                size_t Count = counts[Char];
                size_t index = Char + Count*len; 
                C2 pr = rand_permute[d][index];
                tup_vec.push_back({pr, Char, i});
                counts[Char]++;
            }
            std::sort(tup_vec.begin(), tup_vec.end());
            Vec<std::pair<size_t,C1>> vec;
            for (size_t t=0; t<t_len; t++)  {
                auto tup = tup_vec[t];
                vec.push_back({std::get<2>(tup), std::get<1>(tup)});
            }
            std::sort(vec.begin(), vec.end());
            for (auto pair : vec ) {
                T[d].push_back(pair.second);
            } 
        }
    }
    void sketch_join(const Vec<C1> &seq, Vec<C2>  &T, int sig_len, int t_len) {
        Vec2D<C2> Tvec;
        sketch(seq, Tvec, sig_len, t_len);
        size_t dim = rand_permute.size();
        size_t uniq_siglen = rand_permute[0].size();
        T = Vec<C2>(dim, 0);
        for (size_t d=0; d<dim; d++ ) {
            for (size_t t=0; t<Tvec[d].size(); t++) {
                T[d] += T[d]*uniq_siglen + Tvec[d][t];
            }
        }
    }


    string_opts ops;
    omp_sketch(const string_opts ops) : ops(ops) {}

    void init_permute() {
        init_permute(ops.sig_len, ops.len, ops.dim) ;
    }

    void sketch(const Vec<C1> &seq, Vec2D<C2> &T) {
        sketch(seq, T, ops.sig_len, ops.t_len);
    }
    void sketch_join(const Vec<C1> &seq, Vec<C2> &T) {
        sketch_join(seq, T, ops.sig_len, ops.t_len);
    }
};
