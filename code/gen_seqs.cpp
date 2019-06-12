#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
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
    int gene_len = 100, gene_len2 = 200 , num_genes = 5, 
        padding = 200, num_seq = 5, num_seq2 = 5,repeat = 2;
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
    }
};


void gen_seqs(seq_options opt) {
}

int main(int argc, char* argv[] ) 
{
    seq_options opt; 
    opt.read_args(argc, argv);

    cout << "opt.mutation_rate = " <<  opt.mutation_rate << endl; 
    cout << " save dir = " << opt.save_directory << endl;
    for (int i = 0; i<5; i++) {
        std::srand(std::time(nullptr)); // use current time as seed for random generator
        int random_variable = std::rand();
        cout << "rand int = " << random_variable << endl;
    }

}
