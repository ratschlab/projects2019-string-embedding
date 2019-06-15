#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
using std::cout;
using std::endl;
using std::vector; 
using std::string;


int make_directory(std::string path) {
    cout << " path =  "<< path << endl;
    int out = mkdir(path.c_str(), S_IRWXU);
    return out;
}


#include <cstdlib>
class config {
public:
    std::string file_name,
                fasta_file, val_file, maf_file,
                home_dir,
                project_dir,
                config_path,
                result_path,
                data_path, 
                index_path,
                search_path,
                eval_path,
                log_path;
    void load_proj_dir() {
        home_dir = std::getenv("HOME");
        config_path = home_dir + "/.config_string_embedding";
        std::ifstream fc(config_path);
        if (! fc) {
            std::cerr << " can't open the config file " << endl;
            exit(1);
        }
        string line;
        while (std::getline(fc, line)) {
            cout << " line = " << line << endl;
            if ( line.find("PROJ_DIR") != std::string::npos) {
                project_dir = string(line.begin() + line.find("=") + 1, line.end());
                int c = 0;
                while ( project_dir[c]==' ')
                    c ++;
                project_dir = string(project_dir.begin() + c, project_dir.end());
                cout << " PROJ_DIR = " << project_dir << endl;
                data_path = project_dir + "/data/fasta/" + file_name;
                fasta_file = data_path + "/seqs.fa";
                val_file = data_path + "/seqs.val";
                maf_file = data_path + "/seqs.maf";
                return;
            }
        }
        std::cerr << " couldn't find PROJ_DIR in the config file " << std::endl;
        exit(1);
    }
    config(std::string file_name, std::string result_dir) : file_name(file_name) {
        load_proj_dir();
    }
    config(std::string file_name) : file_name(file_name)  {
        load_proj_dir();
    }
};



struct seq_options {
    bool check_arg(char* str, string name1, string name2) {
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
    string file_name = "tmp";
    

    void read_args(int argc, char* argv[]) {
        using std::stoi;
        using std::stof;
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
            if (check_arg(argv[i],"-f", "--file-name")) {
                file_name = string(argv[++i]);
            } else {
                std::cerr << "illegal argument " << argv[i] << std::endl;
            } 
        }
        // if not provided they default to the other value
        if (num_seq2<0)
            num_seq2 = num_seq;
        if (gene_len2<0)
            gene_len2 = gene_len;
    }
};

string vec2str(vector<char> &vec) {
    string str(vec.begin(), vec.end());
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



void rand_seq(vector<char> &seq, string Sigma, int len) {
    int sigma = Sigma.length();
    for (int i=0; i<len; i++) 
        seq.push_back( Sigma[rand_int(sigma)] );
}

vector<char> rand_seq(string Sigma, int len) {
    vector<char> vec;
    vec.reserve(len);
    rand_seq(vec, Sigma, len);
    return vec;
}


int  mutate_seq(vector<char> &seq, vector<int> &val, vector<char> &seq2, vector<int> &val2, seq_options &opts, string Sigma) {
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
                    val2.push_back(0);
                    i++ ;
                } 
                // insert
                else if (ri==1) {
                    seq2.push_back(Sigma[rand_int(4)]);
                    val2.push_back(0);
                } 
                // delete
                else if (ri==2) { 
                    i++; 
                }
            }
        } else {
            seq2.push_back(seq[i]);
            val2.push_back(val[i]);
            i++;
        }
    }
    return seq2.size() - init_size;
}

/*
vector<char> mutate_seq(vector<char> &seq, vector<int> &val, seq_options &opts, string Sigma) {
    vector<char> seq2;
    mutate_seq(seq, seq2, opts, Sigma);
    return seq2;
}
*/


struct gene_interval {
public:
    int gene_code, begin, end;
    gene_interval(int g, int b, int e ) : gene_code(g), begin(b), end(e) {}
};


void generate_seqs(vector<vector<char> >  &seqs, 
                   vector<vector<int> > &vals, 
                   vector<vector< gene_interval> > &all_intervals, 
                   seq_options &opt) {
    string Sigma = "acgt";
    auto gene_lens = rand_vect(opt.gene_len, opt.gene_len2+1, opt.num_genes);
    auto num_seqs = rand_vect(opt.num_seq, opt.num_seq2+1, opt.repeat);
    
    for (int ri=0; ri< opt.repeat; ri++) {
        vector<vector<char> > genes;
        vector<vector<int> > gene_vals(opt.num_genes, vector<int> ());
        for (int k=0; k<opt.num_genes; k++) {
            genes.push_back(rand_seq(Sigma, gene_lens[k]));
            gene_vals[k].resize(gene_lens[k]);
            std::iota(gene_vals[k].begin(), gene_vals[k].end(), 1);
        }
        cout << " gene_vals = " << endl;
        for (int k=0; k<opt.num_genes; k++) {
            for (int l=0; l<gene_lens[k]; l++) {
                cout << gene_vals[k][l];
            }
            cout << endl;
        }
        cout << " ******* " << endl;

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
                    val.insert(val.end(), gene_vals[k].begin(), gene_vals[k].end());
                } else {
                    gene_len = mutate_seq(genes[k], gene_vals[k],  seq, val,  opt, Sigma); 
                }
                int end = seq.size();
                gene_interval_vec.push_back(gene_interval(gene_code, begin, end));
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
    config conf(opt.file_name);
    make_directory(conf.data_path);
    std::ofstream ffasta, fvals, fmaf;
    ffasta.open (conf.fasta_file);
    fmaf.open (conf.maf_file);
    ffasta << " test " << endl;
    fmaf << " test " << endl;
    //fvals.open (conf.val_file);

    for (int i=0; i<seqs.size() ; i++) {
        auto str = vec2str(seqs[i]);
        cout << str << endl;
        for (int j=0; j<vals[i].size(); j++)
                cout << vals[i][j];
        cout << endl;
    }
    ffasta.close();
    fmaf.close();

    return 0;
}
