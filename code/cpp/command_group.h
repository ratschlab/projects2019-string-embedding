#ifndef __COMMAND_GROUP__ 
#define __COMMAND_GROUP__ 

#include <iostream>
#include <string>
#include <vector>
#include <map>


struct command_group {
    typedef std::string string;

    struct cmd_type ;
    struct cmd_type_int;
    struct cmd_type_str ;
    struct cmd_type_bool;
    
    std::vector<cmd_type*> G;
    std::map<string, int> S, L;
    int __argc; 
    char ** __argv;

    void add(cmd_type *C) ;

    void print_help() ;

    string get_config(int __long = 0) ;

    cmd_type* read(string key) ;

    virtual void add_options() = 0;
    
    bool check_arg(string str, string name1, string name2) ;

    void add_int(int &val, string sh, string lo, string m) ;
    void add_str(string &val, string sh, string lo, string m) ;
    void add_bool(bool &val, bool set, string sh, string lo, string m) ;
    void read_args(int argc, char* argv[]) ;

    void deep_copy(const command_group &cg ) ;
};

struct string_opts : public command_group {
    int k_len = 4, t_len = 2, len = 100, sig_len = 4, dim = 50, num_bins = 5, disc = 8; 
    bool normalize = false;
    string dir = "tmp", res_path = dir + "/res.txt", conf_path = dir+"/conf.txt", src;

    void add_options();

    string_opts(const string_opts &);

    string_opts();
};

#endif //  __COMMAND_GROUP__ 
