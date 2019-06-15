#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
#include <tuple>
using namespace std;



struct seq_options {
    bool check_arg(char* str, std::string name1, std::string name2) {
        if (name1.compare(str)==0 or name2.compare(str)==0)
            return true;
        else
            return false;
        }
public:
    double mutation_rate = .1, reverse_p = 0, geometric_p = .4;
    int gene_len = 100, gene_len2 = -1 , num_genes = 5, 
        padding = 200, num_seq = 5, num_seq2 = -1,repeat = 2;
    bool colored = false, print_values = false, print_seqs = false;
    string save_directory = "tmp.fa";

    void read_args(int argc, char* argv[]) {
        for (int i=1; i<argc; i++) {
            if (check_arg(argv[i],"-m", "--mutation-rate")) {
                mutation_rate = stof(argv[++i]);
            } else  
            if (check_arg(argv[i],"-g", "--gene-len")) {
                gene_len = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-G", "--num-genes")) {
                num_genes = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-P", "--padding")) {
                padding = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-n", "--num-seq")) {
                num_seq = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-r", "--repeat")) {
                repeat = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-V", "--reverse-p")) {
                reverse_p = stof(argv[++i]);
            } else 
            if (check_arg(argv[i],"-X", "--gene-len2")) {
                gene_len2 = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-Y", "--num-seq2")) {
                num_seq2 = stoi(argv[++i]);
            } else 
            if (check_arg(argv[i],"-Z", "--geometric-p")) {
                geometric_p = stof(argv[++i]);
            } else 
            if (check_arg(argv[i],"-c", "--colored")) {
                colored = true;
            } else 
            if (check_arg(argv[i],"-v", "--print-values")) {
                print_values = true;
            } else 
            if (check_arg(argv[i],"-p", "--print-seqs")) {
                print_seqs = true;
            } else 
            if (check_arg(argv[i],"-d", "--save-directory")) {
                save_directory = string(argv[++i]);
            } else {
                std::cout << "illegal argument " << argv[i] << std::endl;
            } 
        }
        // if not provided they default to the other value
        if (num_seq2<0)
            num_seq2 = num_seq;
        if (gene_len2<0)
            gene_len2 = gene_len;
    }
};

std::string vec2str(std::vector<char> &vec) {
    std::string str(vec.begin(), vec.end());
    return str;
}


int rand_int(int min, int max) {
    int rv = std::rand();
    return min + rv% (max - min); 
}

int rand_int(int max) {
    return rand_int(0, max);
}

vector<int> rand_vect(int min, int max, int size) {
    vector<int> vec(size);
    for (int i=0; i<size; i++)
        vec[i] = rand_int(min, max);
    return vec;
}



void rand_seq(vector<char> &seq,string Sigma, int len) {
    int sigma = Sigma.length();
    for (int i=0; i<len; i++) 
        seq.push_back( Sigma[rand_int(sigma)] );
}

vector<char> rand_seq(std::string Sigma, int len) {
    std::vector<char> vec;
    vec.reserve(len);
    rand_seq(vec, Sigma, len);
    return vec;
}


int  mutate_seq(std::vector<char> &seq, std::vector<char> &seq2, seq_options &opts, std::string Sigma) {
    std::default_random_engine gen;
    std::geometric_distribution<int> db(opts.geometric_p);
    int i = 0;
    int init_size = seq2.size();
    while (i< seq.size()) {
        double rv = (double)rand() / RAND_MAX;
        if (rv< opts.mutation_rate) {
            // choose randomly between in/del/sub
            int ri = rand_int(3);         
            // if we've reached the seq end, only insert is possible
            int rounds = db(gen);
            for (int r = 0; r<rounds; r++) {
                if (i==seq.size()) {
                    ri = 1; 
                } 
                // subsitutde
                if (ri==0) {
                    seq2.push_back(Sigma[rand_int(4)]);
                    i++ ;
                } 
                // insert
                else if (ri==1) {
                    seq2.push_back(Sigma[rand_int(4)]);
                } 
                // delete
                else if (ri==2) { 
                    i++; 
                }
            }
        } else {
            seq2.push_back(seq[i]);
            i++;
        }
    }
    return seq2.size() - init_size;
}

std::vector<char> mutate_seq(std::vector<char> &seq, seq_options &opts, std::string Sigma) {
    vector<char> seq2;
    mutate_seq(seq, seq2, opts, Sigma);
}


struct gene_interval {
public:
    int gene_code, begin, end;
    gene_interval(int g, int b, int e ) : gene_code(g), begin(b), end(e) {}
};


void generate_seqs(vector<vector<char>>  &seqs, vector<vector<int> > &vals, vector<vector< gene_interval> > &all_intervals, seq_options &opt) {
    string Sigma = "acgt";
    auto gene_lens = rand_vect(opt.gene_len, opt.gene_len2+1, opt.num_genes);
    auto num_seqs = rand_vect(opt.num_seq, opt.num_seq2+1, opt.repeat);
    
    for (int ri=0; ri< opt.repeat; ri++) {
        vector<vector<char> > genes;
        for (int k=0; k<opt.num_genes; k++) {
            genes.push_back(rand_seq(Sigma, gene_lens[k]));
        }
        for (int i=0; i<num_seqs[ri]; i++) {
            vector<char> seq;
            vector<int> val;
            vector<gene_interval> gene_interval_vec; 
            for (int k=0; k<opt.num_genes; k++) {
                int gene_code = (ri*opt.num_genes) + k + 1;
                rand_seq(seq, Sigma, opt.padding);
                val.insert(val.end(), opt.padding, 0); 
                int gene_len;
                int begin = seq.size();
                if (i==0) {
                    seq.insert(seq.end(), genes[k].begin(), genes[k].end());
                    gene_len = gene_lens[k];
                } else {
                    gene_len = mutate_seq(genes[k], seq,  opt, Sigma); 
                }
                int end = seq.size();
                gene_interval_vec.push_back(gene_interval(gene_code, begin, end));
                val.insert(val.end(), gene_len, gene_code);
                rand_seq(seq, Sigma, opt.padding);
                val.insert(val.end(), opt.padding, 0); 

            }
            seqs.push_back(seq);
            vals.push_back(val);
            all_intervals.push_back(gene_interval_vec);
        }

    }
}


int main(int argc, char* argv[] ) 
{
    std::srand(std::time(nullptr)); // use current time as seed for random generator
    seq_options opt; 
    opt.read_args(argc, argv);
    vector<vector<gene_interval> > all_intervals;
    vector<vector<char> > seqs;
    vector<vector<int> > vals;
    generate_seqs(seqs, vals, all_intervals, opt);

    for (int i=0; i<seqs.size() ; i++) {
        auto str = vec2str(seqs[i]);
  //      cout << " seq "  << i <<  " : " << endl;
        cout << str << endl;
  //      for (int j=0; j<vals[i].size(); j++)
  //          if (vals[i][j]>0)
 //               cout << vals[i][j];
                //cout << str[i][j];
//        cout << endl;
    }

    return 0;
}
