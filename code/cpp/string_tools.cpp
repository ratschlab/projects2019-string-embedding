#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <assert.h>
#include <map>

template<typename seq_type,typename int_type>

struct string_tools_t {
    typedef std::vector<seq_type> seq_t;
    typedef std::vector<int_type> ivect_t;
    typedef std::vector<size_t> indices_t;
    typedef std::vector<double> dvect_t;
    typedef std::vector<std::vector<int>> mat2d_int;
    typedef std::vector<std::vector<double>> mat2d_double;
    typedef std::string string; 


    std::default_random_engine gen;


    template <class T>
    T L1(const std::vector<T>  &v1, const std::vector<T> &v2) {
        assert(v1.size()==v2.size());
        T diff = 0;
        for (int i=0; i<v1.size(); i++ ) {
            diff += abs(v1[i]-v2[i]);
        } 
        return diff;
    }

    void translate(string alpha, string &str, seq_t &seq) {
        std::map<char,int> c2i;
        std::map<int,char> i2c;
        for (int i=0; i<alpha.size() ; i++) {
            char c = alpha[i];
            c2i[c] = i;
            i2c[i] = c;
        }
        seq.clear();
        seq.reserve(str.length());
        for(int i=0; i<str.length(); i++) {
            char c = str[i];
            seq.push_back(c2i[c]);
        }
    }

    void ins(seq_t &seq, int pos, int c) {
        seq.insert(seq.begin() + pos, c);
    }
    void del(seq_t &seq, int pos) {
        seq.erase(seq.begin()+pos);
    }
    void sub(seq_t &seq, int pos, int c) {
        seq[pos] = c;
    }
    void rand_op(seq_t &seq,int sig_len, dvect_t prob) {
        assert(prob.size()==3);
        std::discrete_distribution<int> distribution(prob.begin(), prob.end());
        std::uniform_int_distribution<int> pos_dist(0,seq.size()-1), char_dist(0,sig_len-1);
        int op = distribution(gen), pos = pos_dist(gen), c = char_dist(gen);
        switch (op) {
            case 0:  // ins   
                ins(seq,pos,c);
                break;
            case 1: // del
                del(seq,pos);
                break;
            case 2: // sub
                sub(seq,pos,c);;
        }
    }
    void rand_ops(seq_t &seq, int sig_len, int num_ops, dvect_t prob ) {
        for (int i=0; i<num_ops; i++ ) {
            rand_op(seq,sig_len, prob);
        }
    }

    void rand_seq(seq_t &seq, int sig_len, int len) {
        seq.clear();
        seq.reserve(len);
        std::uniform_int_distribution<int> rchar(0, sig_len-1);
        for (int i=0 ; i<len; i++) {
            seq.push_back(rchar(gen));
        }
    }

    void gen_pairs(std::vector<seq_t> &seqs1, std::vector<seq_t> &seqs2, int len, int sig_len, int num_pairs ) {
        std::uniform_int_distribution<int> num_ops(0, len);
        seqs1 = std::vector<seq_t>(num_pairs);
        seqs2 = std::vector<seq_t>(num_pairs);
        for (int i=0; i<num_pairs; i++ ) {
            rand_seq(seqs1[i], sig_len, len);
            seqs2[i] = seqs1[i];
            rand_ops(seqs2[i], sig_len, num_ops(gen), {1.0,1.0,1.0});
        }
    }

    unsigned int edit_distance(const seq_t &s1, const seq_t & s2)
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

    void pseq(const seq_t &seq) {
        for (auto s : seq ) { 
            std::cout << s;
        } std::cout << "\n";
    }

    template <class T>
    void pvec(const std::vector<T> &vec) {
        for (auto v : vec ) {
            std::cout << v << " " ;
        } std::cout << "\n";
    }
    template <class T> 
    void pmat(const std::vector<std::vector<T> > &mat) {
        for (auto row : mat ) {
            for (auto e : row ) {
                std::cout << e << ", ";
            }
            std::cout << "\n";
        }
    }

    template<class T>
    T pow(T base, T p) {
        T r = 1;
        for (int i=0; i<p; i++)  
            r = r * base; 
        return r;
    }

    size_t read_tuple(const seq_t &seq, const  indices_t &ind, int base) {
        size_t r = 0;
        for (auto i : ind  ) { 
            r = (r*base + seq[i]);
        }
        return r;
    }
    bool inc_ind_sorted(indices_t &ind, int len ) {
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

    void seq2kmer(seq_t &seq, ivect_t &kseq, int k, int sig_len ) {
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

};

