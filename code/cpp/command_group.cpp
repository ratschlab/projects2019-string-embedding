#include "command_group.h"

struct command_group::cmd_type {
    string S, L, H, sv;
    void init(string S, string L, string H = "no description " ) {
        this->S = S; this->L = L; this->H = H;
    }
    virtual bool getval() = 0;
    virtual void set(const cmd_type *cd) = 0;
    virtual void set(string s) {};
    virtual string sval() = 0; // type and value 
    virtual void set() {};

    virtual string to_string () {
        return  sval() + "\t" +  S + ",\t" + L + "\t" + H + "\n"; 
    }
};

struct command_group::cmd_type_int : public command_group::cmd_type {
    int &val;
    cmd_type_int(int &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
    void set(string sval) { val = std::stoi(sval); }
    void set(const cmd_type *cd) { val = ((cmd_type_int*)cd)->val; }
    string sval() { return std::to_string(val) + "(int)"; }
    bool getval() { return true;} 
};

struct command_group::cmd_type_str : public command_group::cmd_type {
    string &val;
    cmd_type_str(string &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
    void set(string sval) { val = sval; }
    void set(const cmd_type *cd) { val = ((cmd_type_str*)cd)->val; }
    string sval() { return val + "(str)"; }
    bool getval() { return true;} 
};

struct command_group::cmd_type_bool : public command_group::cmd_type {
    bool &val; bool to;
    cmd_type_bool(bool &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
    bool getval() { return false;} 
    string sval() { return std::to_string(val) + "(bool)"; }
    void set(const cmd_type *cd) { val = ((cmd_type_bool*)cd)->val; }
    void set() { val = to; }
};




void command_group::add(cmd_type *C) {
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

void command_group::print_help() {
    std::cout <<  "SHORT ARG, \t LONG ARG,\t DESCRIPTION\n";
    for (auto C : G) {
        std::cout << (C->S) << "\t\t" << (C->L) <<"\t" << (C->H) << std::endl;
    }
}

command_group::string command_group::get_config(int __long ) {
    string sconf;
    if (__long == 1) {
        sconf =  "VALUE,\t SHORT ARG \t LONG ARG\t DESCRIPTION\n";
        for (auto c : G) {
            sconf += c->to_string()+ "\n";
        }
    } else if (__long == 0 )    {
        for (auto c : G) {
            sconf += c->S+","+c->L + ": " + c->sval() + ", ";
        }
        sconf += "\n"; 
    }  else if (__long == -1 ) {
        for (int argi = 0 ; argi<__argc; argi++  ) {
            sconf += __argv[argi]  ;
            sconf += " " ;
        }
    }
    return sconf;
}

command_group::cmd_type* command_group::read(string key) {
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


bool command_group::check_arg(string str, string name1, string name2) {
    if (name1.compare(str)==0 or name2.compare(str)==0)
        return true;
    else
        return false;
}

void command_group::add_int(int &val, string sh, string lo, string m) {
    add(new cmd_type_int(val,sh,lo,m));
}
void command_group::add_str(string &val, string sh, string lo, string m) {
    add(new cmd_type_str(val,sh,lo,m));
}
void command_group::add_bool(bool &val, bool set, string sh, string lo, string m) {
    add(new cmd_type_bool(val, sh, lo, m));
}

void command_group::read_args(int argc, char* argv[]) {
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
    // save in case 
    __argc = argc;
    __argv = (char**)argv;
}


void command_group::deep_copy(const command_group &cg ) {
    add_options();
    for (auto cp : cg.G ) {
        string key = cp->S;
        G[S[key]]->set(cp);
    }
}


