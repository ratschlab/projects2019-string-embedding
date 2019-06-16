#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <map>
using std::cout;
using std::endl;
using std::vector; 
using std::string;


int make_directory(std::string path) {
    int out = mkdir(path.c_str(), S_IRWXU);
    return out;
}


class config {
public: std::string file_name,
                maf_file,
                data_command, 
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
            if ( line.find("PROJ_DIR") != std::string::npos) {
                project_dir = string(line.begin() + line.find("=") + 1, line.end());
                int c = 0;
                while ( project_dir[c]==' ')
                    c ++;
                project_dir = string(project_dir.begin() + c, project_dir.end());
                cout << " PROJ_DIR = " << project_dir << endl;
                data_path = project_dir + "/data/fasta/" + file_name;
                data_command = data_path + "/command.sh";
                maf_file = data_path + "/MSA.maf";
                return;
            }
        }
        std::cerr << " couldn't find PROJ_DIR in the config file " << std::endl;
        exit(1);
    }
    string fasta_file(int i) {
        return data_path + "/seqs" + std::to_string(i) + ".fa";
    }
    string val_file(int i) {
        return data_path + "/vals" + std::to_string(i) + ".bin";
    } 

    config(std::string file_name, std::string result_dir) : file_name(file_name) {
        load_proj_dir();
    }
    config(std::string file_name) : file_name(file_name)  {
        load_proj_dir();
    }
};


