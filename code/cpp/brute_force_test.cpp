#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>
#include <set>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <map>
#include <sstream>
#include "config.hpp"
#include <bitset>
#include <typeinfo>     // for the `std::type_info` type
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::stoi;
using std::stof;




struct binary_string {
    // internal typedefs, functions, and members

    typedef uint64_t bstr_t;
    typedef std::bitset<16> bs16;
    typedef std::bitset<24> bs24;
    typedef std::bitset<32> bs32;
    typedef std::bitset<48> bs48;
    typedef std::bitset<64> bs64;
    int max(int a, int b) { return (a > b)? a : b;    }
    int min(int a, int b) { return (a > b)? b : a;    }
    int abs_val(int a) { return (a >= 0)? a : (-a);    }
    int step(int a) { return (a >= 0)? a : (0);    }
    vector<vector<int> > inds;

/*
    struct cmd {
        string S, L, H;
        void init(string S, string L, string H = "no description " ) {
            this->S = S; this->L = L; this->H = H;
        }
        virtual void set(string s)  = 0; 
    };

    struct cmdint : public cmd {
        int &val;
        cmdint(int &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
        void set(string sval) { val = stoi(sval); }
    };
    
    struct cmdstr : public cmd {
        string &val;
        cmdstr(string &v, string S, string L, string H ) : val (v){ this->init(S,L,H); }
        void set(string sval) { val = sval; }
    };

    struct cmd_help : public cmdstr {
        cmd_help(string &s) :  
            cmdstr (s,"-h", "--help", " to get the list of options") { }
    };

    struct opts_t {
        string emp;
        vector<cmd*> G;
        map<string, int> S, L;
        cmd_help *ch = new cmd_help(emp);

        opts_t () {
            add(ch);
        }
        
        void add(cmd *C) {
            string Sh = C->S, Lo = C->L;
            if ( S.count(Sh)>0) {
                cerr << "option (" << Sh << ") already exists, not added " << endl; 
            } else if (L.count(Lo)>0) {
                cerr << "option (" << Lo << ") already exists, not added " << endl; 
            }
            S[C->S] = G.size(); 
            L[C->L] = G.size();
            G.push_back(C);
        }

        void help() {
            cout << "val\t short arg \t long arg\t description\n";
            for (auto C : G) {
                cout << (C->S) << "\t" << (C->L) <<"\t" << (C->H) << endl;
            }
        }

        void read(string key,string sval) {
            cmd* Cp;
            
            if ( S.count(key)==0 and L.count(key)==0) {
                cerr << "option (" << key << ") is not valid " << endl; 
            } else if (S.count(key)>0) { 
                Cp = G[S[key]];  
            } else { 
                Cp = G[L[key]]; 
            } 
            if (Cp==ch) {
                help();
                exit(0);
            }
            //cout << Cp << endl;
            Cp->set(sval);
        }

    };
    */

public:
    int len = 10, sig_len = 2, k_len = 2, t_len = 2, nbits =1, sth = 0;
    int kbits, tbits, ksize, tsize, num;

    void print_params() {
    }

    void update() {
        kbits = nbits * k_len;
        tbits = kbits * t_len;
        ksize = (1<<kbits);
        tsize = (1<<tbits);
        num = (1<<(len*nbits));
    }

    bool check_arg(char* str, string name1, string name2) {
        if (name1.compare(str)==0 or name2.compare(str)==0)
            return true;
        else
            return false;
    }
    
    void read_args(int argc, char* argv[]) {
        /*
        opts_t ops;
        ops.add(new cmdint(k_len,"-k","--kmer-size","size of kmer"));
        ops.add(new cmdint(nbits,"-B","--num-bits","number of bits per character"));
        ops.add(new cmdint(sth,"-S","--sim-thresh","threshold for similarity of kmers"));
        ops.add(new cmdint(t_len,"-t","--tuple-size","number of elemtns in the tuple"));
        ops.add(new cmdint(len,"-N","--str-len","lenght of string to search"));
        ops.add(new cmdint(sig_len,"-s","--sigma-len","size of the alphabet"));
        for (int i=1; i<argc; i++) {
            string key(argv[i]), val(argv[++i]);
            ops.read(key, val);
        }
        exit(0);
       */ 
        
        for (int i=1; i<argc; i++) {
            if (check_arg(argv[i],"-k", "--kmer-size")) {
                k_len = stof(argv[++i]);
            } else if (check_arg(argv[i],"-B", "--num-bits")) {
                nbits = stoi(argv[++i]);
            } else if (check_arg(argv[i],"-S","--sim-threshold")) {
                sth = stoi(argv[++i]);
            } else if (check_arg(argv[i],"-t", "--tuple-size")) {
                t_len = stoi(argv[++i]);
            } else if (check_arg(argv[i],"-N", "--str-len")) {
                len = stoi(argv[++i]);
            } else if (check_arg(argv[i],"-s", "--sig-len")) {
                sig_len = stoi(argv[++i]);
            } else {
                std::cerr << "illegal argument " << argv[i] << std::endl;
            }
        }
        update();
    }



    bstr_t digit(bstr_t &X, int i) {
        bstr_t mask = ((1<<nbits)-1);
        bstr_t d = (X >> (nbits*i));
        d = (d & mask); 
        return d;
    }

    bstr_t kmer(bstr_t &X, int i) {
        bstr_t mask = ((1<<(nbits*k_len))-1);
        bstr_t d = (X >> (nbits*i));
        d = (d & mask); 
        return d;
    }

    bstr_t kmers(bstr_t &X, vector<int> inds) {
        bstr_t ds = 0;
        for (auto i : inds ) {
            ds = (ds<<(nbits*k_len)) + kmer(X,i);
        }
        return ds;
    }

    bstr_t ins(bstr_t X, int i, int c) {
        int bi = nbits*i;
        bstr_t mask = ((1<<nbits)-1);
        bstr_t Xh = ((X>>bi)<<bi), Xl = X-Xh;
        X = (Xh<<nbits) + (c<<bi) + Xl;
        return X;
    }

    bstr_t sub(bstr_t  X, int i, int c) {
        bstr_t mask = ~(((1<<nbits)-1) << (nbits*i));
        X = (X & mask) + (c<<(nbits*i));
        return X;
    }



    int lcs( bstr_t X, bstr_t Y, int m , int n ) {
        vector<int> L1(n+1,0);
        for (int i = 0; i <= m; i++) {
            vector<int>  L2(n+1,0);
            for (int j = 0; j <= n; j++) {
                bstr_t x = digit(X,i-1), y = digit(Y,j-1);
                if (i == 0 || j == 0) {
                    L2[j] = 0;
                } else if (x == y) {
                    L2[j] = L1[j-1] + 1;
                } else {
                    L2[j] = max(L1[j],L2[j-1]);
                }
            }
            L1.swap(L2);
        }
        return L1[n];
    }

    int lcs(bstr_t X, bstr_t Y) {
        return lcs(X,Y,len,len);
    }

    void sorted_inds(int n) {
        vector<int> cur(t_len+1, 0);
        while (cur[t_len]==0) {
            int i=0;
            cur[i]++;
            while (cur[i]>=(n-k_len+1)) {
                cur[i] = 0;
                cur[i+1] ++;
                i++;
            }
            bool sorted = true;
            for (int j=0; j<t_len-1; j++)
                if (cur[j]>=cur[j+1]) {
                    sorted=false; 
                    break;
                }
            if (sorted) {
                inds.push_back(vector<int>(cur.begin(), cur.end()-1));
            }

        }
    }

    void tensor(bstr_t X, vector<int> &T) {
        for (auto ind : inds) {
            bstr_t digs = kmers(X, ind);
            if (digs>=tsize) {
                std::cout << "digs, tsize = " << digs << " " << tsize << endl;
                std::cerr << "digs out of tuple range " << endl;
            }
            T[digs] ++;
        }
    }


    int tensor_diff(vector<int> &Tx, vector<int> &Ty) {
        if (Tx.size() != Ty.size() ) {
            std::cerr << "ERROR: Tx.size() != Ty.size()" << endl;
        }
        int diff = 0;
        for (int i = 0; i<tsize; i++) {
            diff += abs_val(Tx[i]-Ty[i]);
        }
        return diff;
    }

    void simulate () {
        vector<int> Dmax(len+1,-1), Dmin(len+1,1000);
        vector<int> Xm = Dmin, Ym = Dmin;
        vector<vector<int> > Ts(num,vector<int>(tsize,0));
        vector<vector<int> > TW(num,vector<int>(tsize,0));
        sorted_inds(len-k_len+1);
        for (bstr_t X=0; X<num; X++) {
            tensor(X,Ts[X]);
        }
        // compute kmer similarity
        /*
        vector<vector<int>> kmer_sim(ksize,vector<int>(ksize));
        for (int i=0; i<ksize; i++) {
            for (int j=0; j<ksize; j++) {
                kmer_sim[i][j] = lcs(i,j,k_len, k_len) ;
            }
        }
        bstr_t mask = kbits-1;
        vector<vector<int>> tup_sim(tsize,vector<int>(tsize,1));
        for (int ti=0; ti<tsize; ti++) {
            for (int tj=0; tj<tsize; tj++) {
                int s = 0;
                for (int tn=0; tn<t_len; tn++) {
                    bstr_t X = (ti>>(kbits*tn)), Y = (tj>>(kbits*tn));
                    X &= mask; 
                    Y &= mask;
                    s += kmer_sim[X][Y];
                }
                tup_sim[ti][tj] = step(s- sth);
            }
        }

        for (bstr_t X=0; X<num; X++) {
            for (int ti=0; ti<tsize; ti++) {
                for (int tj=0; tj<tsize; tj++) {
                    TW[X][ti] += tup_sim[ti][tj]*Ts[X][tj];
                }
            }
        }
        */
        for (bstr_t X=0; X<num; X++) {
            for (bstr_t Y=X+1; Y<num; Y++) {
//                int l = lcs(X,Y), D = tensor_diff(TW[X],TW[Y]);
                int l = lcs(X,Y), D = tensor_diff(Ts[X],Ts[Y]);
                if (D<Dmin[l]) {
                    Xm[l] = X;
                    Ym[l] = Y;
                }
                Dmin[l] = min(Dmin[l],D);
                Dmax[l] = max(Dmax[l],D);
            }
        }
        Dmin[len] = 0; Dmax[len] = 0;

        cout << "l\t Dmin\t Dmax\tX\t Y\n";
        for (int l=1; l<=len; l++) {
            cout << l << "\t" << Dmin[l] << "\t" << Dmax[l] << 
                "\t" << bs16(Xm[l]) <<"\t" << bs16(Ym[l]) << endl;
        }

    }
};


int main(int argc, char* argv[]) {
    binary_string bst;
    bst.read_args(argc,argv);
    bst.simulate();

    return 0;
}
