#include <string>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <random>
#include <sstream>
#include <cmath>
#include <fstream>
using namespace std;

void print_vec(vector<long int> &vec, string name) {
    cout << name << " = \t";
    for (long int v : vec ) {
        cout << v << ", ";
    }
    cout << endl;
}

void print_svec(vector<long int> &svec, string name, string Sigma) {
    cout << name << " = \t";
    for (long int v : svec ) {
        if (v>=Sigma.size()) {
            cerr << "svec value " << v << " is out of bound for language size = " << Sigma.length() << endl;
            exit(-1);
        }
        cout << Sigma[v];
    }
    cout << endl;
}

struct kmer_tools { 
    string Sigma;
    long int sigma_len, k, Ksigma_len;
    vector<long int> c2i;
    vector<long int> counter;

    kmer_tools() {}

    kmer_tools(string Sigma, long int k) : Sigma(Sigma), sigma_len(Sigma.length()), k(k) {
        c2i = vector<long int>(256, -1);
        for (long int i=0; i<sigma_len; i++) {
            c2i[(long int)Sigma[i]] = i;
        }
        Ksigma_len = sigma_len;
        for (long int i=1; i<k; i++)
            Ksigma_len *= sigma_len;
    }
    long int cindex(char c) {
        long int i =  c2i[(uint8_t)c];
        if (i<0) {
            cerr << "character " << c << " not in the alphabet: " << Sigma << endl;
            exit(-1);
        }
        return i;
    }

    void str2svec(vector<long int> &svec, string &s) {
        svec.clear();
        svec.reserve(s.length());
        for (long int i=0; i<s.length(); i++)
            svec.push_back(s[i]);
    }

    string svec2str(vector<long int> &svec) {
        string s = "";
        for (auto v : svec) {
            s = s + Sigma[v];
        }
        return s;
    }

    void svec2kmers(vector<long int> &kmers,vector<long int> &s) {
        long int slen = s.size();
        kmers.clear();
        kmers.resize(slen-k+1);
        for (long int i=0; i<slen-k+1; i++) {
            kmers[i] = 0;
            for (long int j=0; j<k; j++) {
                kmers[i] = kmers[i]*sigma_len + s[i+j];
            }
        }
    }

    void str2kmers(vector<long int> &kmers,string &s) {
        long int slen = s.length();
        kmers.clear();
        kmers.resize(slen-k+1);
        for (long int i=0; i<slen-k+1; i++) {
            kmers[i] = 0;
            for (long int j=0; j<k; j++) {
                kmers[i] = kmers[i]*sigma_len + cindex(s[i+j]);
            }
        }
    }

    void unify_kmers(vector<long int> &ukmers, vector<long int> &kmers) {
        long int len = kmers.size();
        counter.clear();
        counter.resize(Ksigma_len);
        std::fill(counter.begin(), counter.end(), 0);
        ukmers.clear();
        ukmers.resize(len);

        for (long int i=0; i<len; i++) {
            ukmers[i] = counter[kmers[i]] * Ksigma_len + kmers[i]; 
            counter[kmers[i]] ++;
        }
    }

    void vec2pairs(vector<long int> &pairs, vector<long int> ukmers, long int max_len) {
        long int len = ukmers.size();
        pairs.clear();
        pairs.reserve(len*(len-1)/2);
        long int ulen = Ksigma_len * max_len;
        for (long int i=0; i<len; i++ ) { 
            for (long int j=i+1; j<len; j++) {
                pairs.push_back(ukmers[i]*ulen + ukmers[j]);
            }
        }
    }

    long int pair_diff(vector<long int> &p1, vector<long int> &p2) {
        sort(p1.begin(), p1.end());
        sort(p2.begin(), p2.end());
        long int l1 = p1.size(),
            l2 = p2.size(), 
            i= 0 , j =0,c = 0; 
        while (i<l1 and j<l2 ) {
            if (p1[i]==p2[j]) {
                c++;
                i++;
                j++;
            } else if (p1[i] < p2[j]) {
                i++;
            } else {
                j++;
            }
        }
        long int d = l1 + l2 - 2*c;
        d += (l1-i) + (l2-j);
        return d;
    }
};

struct string_distance {
    string Sigma;
    kmer_tools ktools;
    vector<long int> kmers, ukmers, pairs, svec, proj;
    vector<long int> kmers2, ukmers2, pairs2, svec2, proj2;
    long int k ;
    bool verbose;

    long int len_max, dim;
    vector<vector<long int> > rmat;

    string_distance(string Sigma, long int k, bool verbose = false) : Sigma(Sigma), k(k), verbose(verbose) {
        ktools = kmer_tools(Sigma, k);
    }

    string svec2str(vector<long int> &svec) {
        return ktools.svec2str(svec);
    }

    void set_rmat(long int dim, long int lmax) {
        len_max = lmax*ktools.Ksigma_len; 
        len_max = len_max * len_max; // to maake a pair
        cout << " len_max = " << len_max << endl;
        this-> dim = dim;

        rmat = vector<vector<long int>>(len_max, vector<long int>(dim, 0));
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dis(0, 1);

        cout << " setting rmat " << endl;
        for (long int i=0; i<len_max; i++)
            for (long int j=0; j<dim; j++) 
                rmat[i][j] = 2*dis(gen) - 1;
    }

    void rproj_pair(vector<long int> &pair, vector<long int> &proj) {
        proj.resize(dim, 0);
        for (long int d = 0; d<dim; d++) {
            for (long int i=0; i<pair.size(); i++) {
                proj[d] += rmat[pair[i]][d];
            }
        }
    }

    long int l1(vector<long int> &v, vector<long int> &v2) {
        long int d = 0;
        for (long int i=0; i<v.size(); i++)
            d += abs(v[i]-v2[i]);
        return d;
    }

    long int l2(vector<long int> &v, vector<long int> &v2) {
        long long int d = 0, dl;
        for (long int i=0; i<v.size(); i++) {
            dl =(v[i]-v2[i]);
            d += dl*dl;
            if (d<0) {
                cout << " d before " << (d-dl*dl) << " d now " << d << endl;
                exit(1);
            }
        }
        return d;
    }



    long int inversions(vector<long int> &s1, vector<long int> &s2, long int &L1, long int &L2) {
        long int max_len = max(s1.size(), s2.size());
        ktools.svec2kmers(kmers, s1); 
        ktools.unify_kmers(ukmers, kmers); 
        ktools.vec2pairs(pairs,ukmers, max_len);
        rproj_pair(pairs, proj);
        ktools.svec2kmers(kmers2, s2); 
        ktools.unify_kmers(ukmers2, kmers2); 
        ktools.vec2pairs(pairs2,ukmers2, max_len);
        rproj_pair(pairs2, proj2);
        long int d = ktools.pair_diff(pairs, pairs2);
        L1 = l1(proj, proj2);
        L2 = l2(proj, proj2);
        if (verbose) {
            print_svec(s1, "s1", Sigma);
            print_svec(s2, "s2", Sigma);
            print_vec(kmers, "kmers");
            print_vec(kmers2, "kmers2");
            print_vec(ukmers, "ukmers");
            print_vec(ukmers2, "ukmers2");
            print_vec(pairs, "pairs");
            print_vec(pairs2, "pairs2");
            cout << " pairs diff = " << d << endl;
        }
        return d;
    }

    long int inversions(string &s1, string &s2) {
        ktools.str2svec(svec, s1);
        ktools.str2svec(svec2, s2);
        long int dummy1, dummy2;
        return inversions(svec, svec2, dummy1, dummy2);
    }

    long int min3(long int x, long int y, long int z)
    {
        return min(min(x, y), z);
    }


    long int ed(vector<long int> & str1,vector<long int> &str2)
    {
        long int m = str1.size(), n = str2.size();
        long int dp[m+1][n+1];

        for (long int i=0; i<=m; i++)
        {
            for (long int j=0; j<=n; j++)
            {
                if (i==0)
                    dp[i][j] = j;  // Min. operations = j
                else if (j==0)
                    dp[i][j] = i; // Min. operations = i
                else if (str1[i-1] == str2[j-1])
                    dp[i][j] = dp[i-1][j-1];
                else
                    dp[i][j] = 1 + min3(dp[i][j-1],  // Insert
                            dp[i-1][j],  // Remove
                            dp[i-1][j-1]); // Replace
            }
        }

        return dp[m][n];
    }

    long int ed(string &str1, string &str2)
    {
        long int m = str1.size(), n = str2.size();
        long int dp[m+1][n+1];

        for (long int i=0; i<=m; i++)
        {
            for (long int j=0; j<=n; j++)
            {
                if (i==0)
                    dp[i][j] = j;  // Min. operations = j
                else if (j==0)
                    dp[i][j] = i; // Min. operations = i
                else if (str1[i-1] == str2[j-1])
                    dp[i][j] = dp[i-1][j-1];
                else
                    dp[i][j] = 1 + min3(dp[i][j-1],  // Insert
                            dp[i-1][j],  // Remove
                            dp[i-1][j-1]); // Replace
            }
        }

        return dp[m][n];
    }
};

struct str_tools {
    string Sigma;
    long int sigma_len;
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()

    str_tools(string Sigma) : Sigma(Sigma), sigma_len(Sigma.length()), gen(std::mt19937(rd())){}
    str_tools(string Sigma, long int seed) : Sigma(Sigma), sigma_len(Sigma.length()), gen(std::mt19937(seed)){}


    void print(vector<long int> &svec) {
        for (long int v : svec)
            cout << Sigma[v];
        cout << endl;
    }

    void rand_svec(vector<long int> &svec, long int len) {
        std::uniform_int_distribution<> dis(0, sigma_len-1);

        svec.clear();
        svec.resize(len);
        for (long int i=0; i<len; i++) {
            svec[i] = dis(gen); 
        }
    }

    void rand_ed_one(vector<long int> &svec) {
        long int slen = svec.size();
        std::uniform_int_distribution<> dis(0, 2);
        long int choice = dis(gen);
        if (choice==0) { // insert 
            std::uniform_int_distribution<> dis(0, sigma_len-1), dis2(0, slen);
            long int nc = dis(gen), pos = dis2(gen); 
            svec.insert(svec.begin()+pos, nc);
        } else if (choice ==1) { // delete
            std::uniform_int_distribution<> dis(0, sigma_len-1), dis2(0, slen-1);
            long int pos = dis2(gen); 
            svec.erase(svec.begin()+pos, svec.begin()+pos+1);
        } else { // substitude
            std::uniform_int_distribution<> dis(0, sigma_len-1), dis2(0, slen-1);
            long int nc = dis(gen), pos = dis2(gen); 
            while (nc==svec[pos])
                nc = dis(gen);
            svec[pos] = nc;
        }
    }

    void edit(vector<long int> &svec, long int ed) {
        for (long int i=0; i<ed; i++)
            rand_ed_one(svec);
    }

    void make_rand_pair(vector<long int> &svec, vector<long int> &svec2, long int len, long int ed) {
        rand_svec(svec, len);
        svec2 = svec;
        edit(svec2, ed);
    }

};



int main ( int argc, char* argv[]) {
    long int N = 100, k = 4, len = 100, dim = 20;
    string fname = "tmp";
    for (long int i=1; i<argc; i++ ) {
        if (strcmp(argv[i],"-k")==0) {
            k = stoi(argv[++i]);
        } else if (strcmp(argv[i],"-l")==0) {
            len = stoi(argv[++i]);
        } else if (strcmp(argv[i],"-N")==0) {
            N = stoi(argv[++i]);
        } else if (strcmp(argv[i],"-o")==0) {
            fname = argv[++i];
        } else if (strcmp(argv[i],"-d")==0) {
            dim = stoi(argv[++i]);
        }
    }
    fname = "results/" + fname;
    fname = fname + "_N" + to_string(N) + "_k" + to_string(k) + "_L" + to_string(len) + "_D" + to_string(dim)+ ".txt";
    // init
    string Sigma = "acgt";
    bool verbose = false;
    string_distance sdist(Sigma, k, verbose);
    sdist.set_rmat(dim,len); 
    str_tools st(Sigma);
    vector<long int> svec, svec2;
    ofstream fout;
    fout.open (fname);
    long int L1, L2;
    for (long int ed = 1; ed<=len; ed++) {
        cout << " ed " << ed << endl;
        for (long int i=0; i<N; i++) {
            st.make_rand_pair(svec, svec2, len, ed);
            long int ed2 = sdist.ed(svec,svec2); 
            long int ic = sdist.inversions(svec,svec2, L1, L2);
            fout << ed2 << ", " << ic << ", " << L1 << ", " << L2 << endl;
//            cout << L2 << " ";
        }
    }
    return 0;
}
