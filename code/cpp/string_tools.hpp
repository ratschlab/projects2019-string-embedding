#ifndef __STRING_TOOLS__
#define  __STRING_TOOLS__


#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <cassert>
#include <map>
#include "string_options.hpp"
#include "st_utils.hpp"

struct string_tools_t : public st_utils {
    typedef std::string string;
    template<typename T>
    using Vec = std::vector<T>;
    template<typename T>
    using Vec2D = Vec<Vec<T>>;
    template<typename T>
    using Vec3D = Vec2D<Vec<T>>;
    std::default_random_engine gen;


    template<typename T>
    void translate(const string &alpha, const string &str, Vec<T> &seq) {
        std::map<char, T> c2i;
        std::map<T, char> i2c;
        for (int i = 0; i < alpha.size(); i++) {
            char c = alpha[i];
            c2i[c] = i;
            i2c[i] = c;
        }
        seq.clear();
        seq.reserve(str.length());
        for (char c : str) {
            seq.push_back(c2i[c]);
        }
    }

    template<typename T>
    void ins(Vec<T> &seq, int pos, int c) {
        seq.insert(seq.begin() + pos, c);
    }

    template<typename T>
    void del(Vec<T> &seq, int pos) {
        seq.erase(seq.begin() + pos);
    }
    template <typename T> 
        void sub(Vec<T>  &seq, int pos, int c) {
            seq[pos] = c;
        }
    template <typename T> 
        void rand_op(Vec<T> &seq,int sig_len, Vec<double> prob) {
            assert(prob.size()==3);
            std::discrete_distribution<int> distribution(prob.begin(), prob.end());
            std::uniform_int_distribution<int> pos_dist(0,seq.size()-1), char_dist(0,sig_len-1);
            int op = distribution(gen), pos = pos_dist(gen), c = char_dist(gen);
            switch (op) {
                case 0:  // ins   
                    ins(seq, pos, c);
                    break;
                case 1: // del
                    del(seq, pos);
                    break;
                case 2: // sub
                    sub(seq, pos, c);;
                default:
                    exit(1);
            }
        }

    template <typename T> 
        void rand_ops(Vec<T> &seq, int sig_len, int num_ops, Vec<double> prob ) {
            for (int i=0; i<num_ops; i++ ) {
                rand_op(seq,sig_len, prob);
            }
        }

    template <typename T> 
        void rand_seq(Vec<T> &seq, int sig_len, int len) {
            seq.clear();
            seq.reserve(len);
            std::uniform_int_distribution<int> rchar(0, sig_len-1);
            for (int i=0 ; i<len; i++) {
                seq.push_back(rchar(gen));
            }
        }

    template <typename T> 
        void gen_pairs(Vec2D<T> &seqs1, Vec2D<T> &seqs2, int len, int sig_len, int num_pairs ) {
            std::uniform_int_distribution<int> num_ops(0, len);
            seqs1 = Vec2D<T>(num_pairs);
            seqs2 = Vec2D<T>(num_pairs);
            for (int i=0; i<num_pairs; i++ ) {
                rand_seq(seqs1[i], sig_len, len);
                seqs2[i] = seqs1[i];
                rand_ops(seqs2[i], sig_len, num_ops(gen), {1.0,1.0,1.0});
            }
        }

    template <typename T> 
        unsigned int edit_distance(const Vec<T> &s1, const Vec<T> & s2)
        {
            const std::size_t len1 = s1.size(), len2 = s2.size();
            std::vector<std::vector<unsigned int>> d(len1 + 1, std::vector<unsigned int>(len2 + 1));

            d[0][0] = 0;
            for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
            for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

            for(unsigned int i = 1; i <= len1; ++i)
                for(unsigned int j = 1; j <= len2; ++j)
                    d[i][j] = std::min({d[i - 1][j] + 1, d[i][j - 1] + 1, 
                            d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1) });
            return d[len1][len2];
        }


    template <typename T> 
        size_t read_tuple(const Vec<T> &seq, const  Vec<size_t> &ind, int base) {
            size_t r = 0;
            for (auto i : ind  ) { 
                r = (r*base + seq[i]);
            }
            return r;
        }

    bool inc_ind_sorted(Vec<size_t> &ind, int len ) {
        auto cur = ind.size()-1;
        while (ind[cur]  == len-1 - (ind.size()-1 - cur) ) {
            if (cur==0) {
                return false;
            } 
            cur--;
        }
        ind[cur]++;
        for (auto c=cur+1; c<ind.size(); c++) {
            ind[c] = ind[c-1]+1;
        }
        return true;
    }

    template <typename T1, typename T2> 
        void seq2kmer(Vec<T1> &seq, Vec<T2> &kseq, int k, int sig_len ) {
            kseq.clear();
            kseq.reserve(seq.size());
            for (int i=0; i<seq.size()-k+1; i++) {
                int kid = 0;
                for (int j=0; j<k; j++) {
                    kid = kid*sig_len + seq[i+j];
                }
                kseq.push_back(kid);
            }
        }

    template <typename T> 
        void gen_pairs(Vec<T> &seqs1, Vec<T> &seqs2, int num_pairs, const string_opts &ops) {
            gen_pairs(seqs1,seqs2,ops.len, ops.sig_len, num_pairs); 
        }

    template <typename T1, typename T2> 
        void seq2kmer(Vec<T1> &seq, Vec<T2>  &kseq, const string_opts &ops) { 
            seq2kmer(seq, kseq, ops.k_len, ops.sig_len);
        }


};



#endif // __STRING_TOOLS__
