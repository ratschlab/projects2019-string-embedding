#include <iostream>
#include <string>
#include <vector>
#include <map>


struct command_group {
    // command classes 
    typedef std::string string;

    struct cmd_type {
        string S, L, H, sv;
        void init(string S, string L, string H = "no description " ) {
            this->S = S; this->L = L; this->H = H;
        }
        virtual bool getval() = 0;
        virtual void set(string s) {};
        virtual void set(const cmd_type *cd) = 0;
        virtual void set() {}
    };

    struct cmd_type_int : public cmd_type {
        int &val;
        cmd_type_int(int &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
        void set(string sval) { val = std::stoi(sval); }
        void set(const cmd_type *cd) { val = ((cmd_type_int*)cd)->val; }
        bool getval() { return true;} 
    };
    
    struct cmd_type_str : public cmd_type {
        string &val;
        cmd_type_str(string &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
        void set(string sval) { val = sval; }
        void set(const cmd_type *cd) { val = ((cmd_type_str*)cd)->val; }
        bool getval() { return true;} 
    };
    
    struct cmd_type_bool : public cmd_type {
        bool &val; bool to;
        cmd_type_bool(bool &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
        bool getval() { return false;} 
        void set(const cmd_type *cd) { val = ((cmd_type_bool*)cd)->val; }
        void set() { val = to; }
    };
    
    // list of internal memers and functions
    std::vector<cmd_type*> G;
    std::map<string, int> S, L;
    int __argc; 
    char ** __argv;


        
    void add(cmd_type *C) {
        string Sh = C->S, Lo = C->L;
        if ( S.count(Sh)>0) {
            std::cerr << "option (" << Sh << ") already exists, not added " << std::endl; 
        } else if (L.count(Lo)>0) {
            std::cerr << "option (" << Lo << ") already exists, not added " << std::endl; 
        }
        S[C->S] = G.size(); 
        L[C->L] = G.size();
        G.push_back(C);
    }

    void print_help() {
        std::cout << "val\t short arg \t long arg\t description\n";
        for (auto C : G) {
            std::cout << "\t"<< (C->S) << "\t\t" << (C->L) <<"\t" << (C->H) << std::endl;
        }
    }

    cmd_type* read(string key) {
        cmd_type* Cp;
        
        if ( S.count(key)==0 and L.count(key)==0) {
            std::cerr << "option (" << key << ") is not a valid argument " << std::endl; 
            return 0;
        } else if (S.count(key)>0) { 
            Cp = G[S[key]];  
        } else { 
            Cp = G[L[key]]; 
        } 
        return Cp;
    }

    virtual void add_options() = 0;
    
    bool check_arg(string str, string name1, string name2) {
        if (name1.compare(str)==0 or name2.compare(str)==0)
            return true;
        else
            return false;
    }

    void add_int(int &val, string sh, string lo, string m) {
        add(new cmd_type_int(val,sh,lo,m));
    }
    void add_str(string &val, string sh, string lo, string m) {
        add(new cmd_type_str(val,sh,lo,m));
    }
    void add_bool(bool &val, bool set, string sh, string lo, string m) {
        add(new cmd_type_bool(val, sh, lo, m));
    }

    void read_args(int argc, char* argv[]) {
        add_options();

        for (int i=1; i<argc; i++) {
            string cmd_name(argv[i]); 
            if (check_arg(cmd_name,"-h","--help")) {
                print_help();
                exit(0);
            }
            cmd_type *Cm = read(cmd_name);
            if (Cm==0) {
                continue;
            } if  (not Cm->getval()) { 
                Cm->set();
            } else if (i==(argc-1) ) {
                std::cerr << "Command (" << cmd_name << ") requires a value " << std::endl;
            } else {
                Cm->set(argv[++i]);
            }
        }
        // save for deep copy default
        __argc = argc;
        __argv = (char**)argv;
    }


    void deep_copy(const command_group &cg ) {
        read_args(cg.__argc, cg.__argv);
        // override the changes
        for (auto cp : cg.G ) {
            string key = cp->S;
            G[S[key]]->set(cp);
        }
    }
};

