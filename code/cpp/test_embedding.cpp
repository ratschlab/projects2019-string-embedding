#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "st_utils.hpp"
#include "string_options.hpp"
#include "tensor_sketch.hpp"
#include "tensor_embed.hpp"
#include "omp_sketch.hpp"

// TODO seprate test unit from the main file
struct test_unit : public string_tools_t {
    typedef Vec<int> seq_t;
    template<class T>
    using Vec = std::vector<T>;
    template<class T>
    using Vec2D = Vec<Vec<T>>;
    typedef std::string string;

    template<class T>
    std::string seq2string(const Vec<T> &seq) {
        std::string str;
        for (auto s : seq) {
            str += (char) (s + 'a');
        }
        return str;
    }

    void test_pairs(string_opts &opts) {
        string_opts org_opts(opts);
        int num = opts.num_exp;
        st_utils sutil;
        std::ofstream fout(opts.res_path), fconf(opts.conf_path);
        // TODO read project dir from config file
        std::string sys_cmd = "mkdir -p " + opts.dir;
        std::system(sys_cmd.data());
        fconf << opts.get_config(1);
        std::cout << opts.get_config(0);

        Vec<seq_t> S1, S2;
        gen_pairs(S1, S2, num, org_opts);
        tensor_embed<int, double> st(opts);
        opts.sig_len = pow(opts.sig_len, opts.k_len);
        tensor_sketch<int, double> TS(opts);
        tensor_embed<int, double> TE(opts);
        opts.len = (int) (opts.len * 1.5);
        omp_sketch<int, int> OS(opts);
        TS.init_rand();
        TE.init_ideal();
        OS.init_permute();


        std::string header;
        if (opts.verbose) {
            header += "seq1,\t seq2,\t";
        }
        header += "ed, \t td, \t hd, \t HD, OMD \n";
        std::cout << "seq1,\t seq2,\t" << header;
        fout << header;
        for (int i = 0; i < num; i++) {
            seq_t r1 = S1[i], r2 = S2[i];
            Vec<int> s1, s2, T1, T2, osk1, osk2;
            Vec<double> h1, h2, H1, H2;

            st.seq2kmer(r1, s1, org_opts);
            st.seq2kmer(r2, s2, org_opts);

            TE.embed_naive(s1, T1);
            TE.embed_naive(s2, T2);
            TE.sketch_ideal(T1, h1);
            TE.sketch_ideal(T2, h2);
            TS.sketch(s1, H1);
            TS.sketch(s2, H2);
            OS.sketch_join(s1, osk1);
            OS.sketch_join(s2, osk2);

            auto tdiff = sutil.L1_diff(T1, T2);
            auto omh_diff = sutil.hamming_diff(osk1, osk2);
            auto hdiff = sutil.median(h1, h2);
            auto Hdiff = sutil.median(H1, H2);
            auto ed = st.edit_distance(r1, r2);

            std::string input_seqs = seq2string(r1) + "," + seq2string(r2) + ",";
            std::string seq_out;
            if (opts.verbose == 1) {
                seq_out += input_seqs;
            }

            std::cout << seq_out << ed << ",\t " << tdiff << ",\t " <<
                      std::setprecision(4) << hdiff << ",\t " << std::setprecision(4) << Hdiff << ",\t " << omh_diff
                      << "\n";
            fout << input_seqs << ed << ",\t " << tdiff << ",\t " <<
                 std::setprecision(4) << hdiff << ",\t " << std::setprecision(4) << Hdiff << ",\t " << omh_diff << "\n";
        }
    }

    /* 
       int k_len = 4, t_len = 2, len = 100, sig_len = 4, dim = 50, num_bins = 5, disc = 8, num_exp = 100; 
       bool normalize = false;
       string dir = "tmp", res_path = dir + "/res.txt", conf_path = dir+"/conf.txt", src;
       */
    string_opts make_conf() {
        string_opts o;
        o.k_len = 1;
        o.t_len = 2;
        o.len = 4;
        o.sig_len = 2;
        o.dim = 20;
        o.num_bins = 5;
        o.disc = 8;
        o.num_exp = 100;
        o.normalize = false;
        return o;
    }

    void test_omp(const string &str1, const string &str2, const string &alpha) {
        auto opts = make_conf();

        int k = 2, sig_len = alpha.length();
        Vec<int> seq1, seq2;
        Vec<int> kmer1, kmer2;
        translate(alpha, str1, seq1);
        translate(alpha, str2, seq2);
        seq2kmer(seq1, kmer1, k, sig_len);
        seq2kmer(seq2, kmer2, k, sig_len);
        opts.len = 8;
        opts.k_len = 1;
        opts.t_len = 2;
        opts.sig_len = 2;
        omp_sketch<int, int> omp_sk(opts);
        omp_sk.init_permute();
        Vec2D<int> osk1, osk2;
        omp_sk.sketch(seq1, osk1);
        omp_sk.sketch(seq2, osk2);
        pseq(seq1, "seq1");
        pseq(seq2, "seq2");
        pseq(kmer1, "kmer1");
        pseq(kmer2, "kmer2");
        pmat(osk1, "osk1", true);
        pmat(osk2, "osk2", true);
    }
};

int main(int argc, char *argv[]) {
    using std::cout;
    string_opts opts;
    test_unit tu;
    opts.read_args(argc, argv);

    tu.test_pairs(opts);
    return 0;
}
