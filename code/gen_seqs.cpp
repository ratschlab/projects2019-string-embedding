#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
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

vector<char> rand_seq(string Sigma, int len) {
    int sigma = Sigma.length();
    vector<char> seq(len);
    for (int i=0; i<len; i++) 
        seq[i] = Sigma[rand_int(sigma)];
    return seq;
}


vector<char> mutate_seq(vector<char> &seq, seq_options &opts, std::string Sigma) {
    std::default_random_engine gen;
    std::geometric_distribution<int> db(opts.geometric_p);
    vector<char> seq2;
    int i = 0;
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
    string s(seq2.begin(), seq2.end());
    return seq2;
}

void print_seq(vector<char> &vec) {
    for (int i; i<vec.size(); i++)
        cout << vec[i];
    cout << endl;
}



void gen_seqs(seq_options opt) {
    string Sigma = "acgt";
    auto gene_lens = rand_vect(opt.gene_len, opt.gene_len2+1, opt.num_genes);
    auto num_seqs = rand_vect(opt.num_seq, opt.num_seq2+1, opt.repeat);
}

int main(int argc, char* argv[] ) 
{
    std::srand(std::time(nullptr)); // use current time as seed for random generator
    seq_options opt; 
    opt.read_args(argc, argv);
    int len = 100;
    string Sigma = "acgt";
    auto seq = rand_seq(Sigma, len);
    auto seq2 = mutate_seq(seq, opt, Sigma);
    string s1(seq.begin(), seq.end());
    string s2(seq2.begin(), seq2.end());
    cout << " seq1 = " << s1 << endl;
    cout << " seq2 = " << s2 << endl;

    /*
    cout << "opt.mutation_rate = " <<  opt.mutation_rate << endl; 
    cout << " save dir = " << opt.save_directory << endl;
    for (int i = 0; i<5; i++) {
        int random_variable = std::rand();
        cout << "rand int = " << random_variable << endl;
        double f = (double)rand() / RAND_MAX;
        cout << "rand double = " << f << endl; 
    }
    */

}
