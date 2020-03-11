#ifndef __ST_UTILS__
#define __ST_UTILS__

#include <iostream> 
#include <string> 

struct st_utils { 

    template <class T>
    void pseq(const std::vector<T> &seq, std::string name = "") {
        if (name.size()>0) {
            std::cout << name << " = " ;
        }
        for (auto s : seq ) { 
            std::cout << s;
        } std::cout << "\n";
    }

    template <class T>
    void pvec(const std::vector<T> &vec, std::string name = "") {
        if (name.size()>0) {
            std::cout << name << " = " ;
        }
        for (auto v : vec ) {
            std::cout << v << " " ;
        } std::cout << "\n";
    }
    template <class T> 
    void pmat(const std::vector<std::vector<T> > &mat, std::string name = "", bool tr = false) {
        if (name.size()>0) {
            std::cout << name << " = \n" ;
        }
        if (tr==false) {
            for (auto row : mat ) {
                for (auto e : row ) {
                    std::cout << e << ", ";
                }
                std::cout << "\n";
            }
        } else { // transpose  
            for (size_t i=0; i< mat[0].size(); i++  ) {
                for (size_t j=0; j<mat.size(); j++) {
                    std::cout << mat[j][i] << ", " ; 
                }
                std::cout << "\n";
            }
        }
    }

    template<class T, class T2>
    T pow(T base, T2 p) {
        T r = 1;
        for (int i=0; i<p; i++)  
            r = r * base; 
        return r;
    }
    template <class T>
    T L1(const std::vector<T>  &v1) {
        T diff = 0;
        for (int i=0; i<v1.size(); i++ ) {
            diff += abs(v1[i]);
        } 
        return diff;
    }
    template <class T>
    T L1_diff(const std::vector<T>  &v1, const std::vector<T> &v2) {
        assert(v1.size()==v2.size());
        T diff = 0;
        for (int i=0; i<v1.size(); i++ ) {
            diff += abs(v1[i] - v2[i]);
        } 
        return diff;
    }
    template <class T> 
    size_t hamming_diff(const std::vector<T> &v1, const std::vector<T> &v2) {
        assert(v1.size()==v2.size());
        size_t d = 0;
        for (size_t i=0; i<v1.size(); i++ ) {
            d += (v1[i]!=v2[i] ) ? 1 : 0 ;
        }
        return d;
    }


    template <class T>
    T abs(T v) { 
        T r =  (v>=0) ? (v) : (-v); 
        return r;
    }


    template <class T>
   T median(std::vector<T>  v1, const std::vector<T> &v2) {
        assert(v1.size()==v2.size());
        for (int i=0; i<v1.size(); i++ ) {
            v1[i] = abs(v1[i]- v2[i]);
        } 
        std::sort(v1.begin(), v1.end());
        size_t med = v1.size()/2;
        return v1[med];
    }
};


#endif //  __ST_UTILS__
