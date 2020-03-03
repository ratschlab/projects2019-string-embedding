#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include "st_utils.hpp"
#include "string_options.hpp"
#include "tensor_sketch.hpp"
#include "tensor_embed.hpp"


struct test_unit { 
        typedef tensor_sketch<int,double> TE;
        typedef tensor_embed<int,double> TEN;
        typedef string_tools_t ST;
        template < class T> 
        using Vec = std::vector<T>;
        typedef Vec<int> seq_t;
        typedef Vec<int> ivect_t;
        typedef Vec<double> dvect_t;
    void test_pairs(std::string alpha, string_opts &opts) {
        string_opts org_opts(opts);
        int num = opts.num_exp;
        st_utils sutil;
        std::ofstream  fout(opts.res_path), fconf(opts.conf_path);
        std::string sys_cmd = "mkdir -p " + opts.dir;
        std::system(sys_cmd.data());
        fconf << opts.get_config(1);
        std::cout << opts.get_config(0) ;

        Vec<seq_t> S1, S2;
        ST ss;
        TE st(opts);
        opts.sig_len = sutil.pow(opts.sig_len ,opts.k_len);
        TE tt(opts);
        TEN ttn(opts);
        
        tt.init_rand();
        ttn.init_ideal();

        st.gen_pairs(S1,S2,num, org_opts);

        auto header = "ed, \t td, \t hd, \t HD, \n";
        std::cout << header;
        fout << header;
        for (int i=0; i<num; i ++ ){
            seq_t r1 = S1[i], r2 = S2[i]; 
            ivect_t s1, s2, T1, T2;
            dvect_t h1, h2, H1, H2;

            st.seq2kmer(r1,s1, org_opts);
            st.seq2kmer(r2,s2, org_opts);

            ttn.embed_naive(s1,T1);
            ttn.embed_naive(s2,T2);
            ttn.sketch_ideal(T1,h1);
            ttn.sketch_ideal(T2,h2);
            tt.sketch(s1,H1);
            tt.sketch(s2,H2);

            auto tdiff = sutil.L1_diff(T1,T2);
            auto hdiff = sutil.median(h1,h2);
            auto Hdiff = sutil.median(H1,H2);
            auto ed = ss.edit_distance(r1,r2);
            std::cout << ed << ",\t " << tdiff << ",\t " <<  
                std::setprecision(4) <<  hdiff<< ",\t " << std::setprecision(4) << Hdiff<< "\n";
            fout << ed << ",\t " << tdiff << ",\t " <<  
                std::setprecision(4) <<  hdiff<< ",\t " << std::setprecision(4) << Hdiff<< "\n";
        }
    }
};


int main(int argc, char* argv[]) {
    using std::cout;
    string_opts  opts;
    test_unit tu;
    opts.read_args(argc,argv);

    std::string alpha = "acgt"; 
    tu.test_pairs(alpha, opts);

    return 0;
}
