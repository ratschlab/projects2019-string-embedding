#include <iostream>
#include <string> 
#include <vector> 
#include <fstream>
#include <unordered_map>
#include <random>
#include <cmath>
using std::vector;
using std::cout;
using std::endl;
using std::string;

typedef uint8_t kindex_t;
typedef uint64_t index_t;
typedef std::vector<uint8_t>  svec_t;
typedef svec_t::iterator svec_t_it;
typedef std::vector<kindex_t> kvec_t;


struct alphabet {
    std::string Sigma;
    size_t sigma_len;
    std::vector<int> c2i = std::vector<int>(256, -1); 
    std::vector<char> i2c;
    vector<string> kmers;

    alphabet(string Sigma) : Sigma(Sigma) {
        sigma_len = Sigma.length();
        i2c.resize(sigma_len);
        for (int i=0; i< sigma_len; i++) {
            char c = Sigma[i];
            c2i[(int)c] = i;
            i2c[i] = Sigma[i];;
        }
    }

    void translate_str2vec(string &str, svec_t &svec) {
        svec.reserve(str.length());
        for (int i=0; i<str.length(); i++)
            svec.push_back(c2i[(int)str[i]]);
    }

    kindex_t kmer_index(svec_t_it begin, svec_t_it end) {
        size_t i = 0;
        for (auto it = begin; it!=end; it++) {
            i = i*sigma_len + *it;
        }
        return i;
    }
    void translate_svec2kvec(svec_t &svec, kvec_t &dvec, int len) {
        dvec.reserve(svec.size());
        for (int i=0; i<svec.size() - len + 1; i++) {
            auto id = kmer_index(svec.begin() + i, svec.begin() + i + len);
            dvec.push_back(id);
        }
    }

    void all_kmers(int len) {
        kmers = vector<string> () ;
        size_t tlen = 1;
        for (int i=0; i< len; i++)
            tlen = tlen * sigma_len;
        for (int i=0; i<tlen; i++) {
            string kmer = "";
            kmer = kmer + i2c[i % sigma_len];
            size_t len_pow = sigma_len;
            for (int j=1; j<len; j++) {
                kmer = kmer + i2c[ (i/len_pow) % sigma_len];
                len_pow = len_pow * sigma_len;
            }
            kmers.push_back(kmer);
        }
    }

};

struct string_embedding {
    vector<vector<int> > randmat;

    string_embedding(int size) {
        cout << " started embeddgin constructyor " << endl;
        std::random_device r;
        std::default_random_engine e1(0);
        std::uniform_int_distribution<int> uniform_dist(1, 6);
        int mean = uniform_dist(e1);
        cout << " mean = " << mean << endl;
    }

};




int main () {
    string path = "/cluster/work/grlab/projects/projecs2019-string-embedding/synthetic/data/fasta/tmp";
    
    string seq_path = path + "/seqs0.fa";
    std::ifstream infile(seq_path);
    string line, seq;
    std::getline(infile, line);
    std::getline(infile, seq);
    alphabet alpha("actg");
    alpha.all_kmers(2);
    for (auto kmer : alpha.kmers) {
        svec_t svec;
        alpha.translate_str2vec(kmer, svec);
        cout << kmer << " -> " << (int)alpha.kmer_index(svec.begin(), svec.end()) << endl;
    }
    string_embedding se(1);
    return 0;
}
