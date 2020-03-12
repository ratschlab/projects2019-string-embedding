#ifndef __STRING_OPTIONS__ 
#define __STRING_OPTIONS__ 

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "command_group.h"

struct string_opts : public command_group {
    int k_len = 4, t_len = 2, len = 100, sig_len = 4, dim = 50, num_bins = 5, disc = 8, num_exp = 100, verbose = 0; 
    bool normalize = false;
    // TODO read project dir from a config file
    string proj_dir = "/Users/amirjoudaki/Codes/projects2019-string-embedding";
    string dir = proj_dir + "/code/cpp/tmp", res_path = dir + "/res.txt", conf_path = dir+"/conf.txt", src;

    void add_options() override;

    string_opts(const string_opts &);

    string_opts();
};

// TODO seprate sequence generation from analysis
void string_opts::add_options() {
    add_int(k_len, "-k", "--kmer-size", "size of kmer");
    add_int(t_len, "-t", "--tuple-size", "number of elemtns in the tuple");
    add_int(len, "-l", "--str-len", "lenght of string to search");
    add_int(sig_len, "-S", "--sigma-len", "size of the alphabet");
    add_int(verbose, "-v", "--verbose", "verbosity, 1: print sequences 0: only results");
    add_int(dim, "-m", "--dim", "dimension of sketching");
    add_int(num_bins, "-B", "--num-bins", "discretization for tensor sketching");
    add_int(num_exp, "-E", "--num-exp", "number of experiments or tests");
    add_bool(normalize, "-N", "--normaliese", "normlize cauchy output");
    add_str(dir, "-d", "--dest-dir", "normlize cauchy output");
    add_str(src, "-s", "--src-dir", "normlize cauchy output");

}

string_opts::string_opts(const string_opts &so) {
    deep_copy(so);
}

string_opts::string_opts() = default;

#endif // __STRING_OPTIONS__ 
