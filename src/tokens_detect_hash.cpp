#include <Rcpp.h>
#include <vector>
#include <algorithm>  
//#include "quanteda.h"

using namespace Rcpp;
using namespace std;
//using namespace quanteda;


// [[Rcpp::export]]
IntegerVector qatd_cpp_detect_hash_vector(IntegerVector tokens_, 
                                          IntegerVector tokens_loc_, 
                                          IntegerVector seq_,
                                          int id){
    
    std::vector<int> tokens = Rcpp::as< std::vector<int> >(tokens_);
    IntegerVector tokens_loc = tokens_loc_;
    std::vector<int> seq = Rcpp::as< std::vector<int> >(seq_);
    int len_seq = seq.size();
    std::vector<int>::iterator it;
    it = std::search(tokens.begin(), tokens.end(), seq.begin(), seq.end());
    //Rcout << std::distance(tokens.begin(), it) << "\n";
    while(it != tokens.end()){
        int loc = std::distance(tokens.begin(), it);
        tokens_loc[loc] = id;
        std::advance(it, seq.size());
        it = std::search(it, tokens.end(), seq.begin(), seq.end());
        //Rcout << loc << "\n";
    }
    return tokens_loc;
    
}

// [[Rcpp::export]]
List qatd_cpp_detect_hash_list(List texts_, 
                               List texts_loc_,
                               IntegerVector seq,
                               int id){
    
    List texts = texts_;
    List texts_loc = texts_loc_;
    int len = texts.size();
    for (int h = 0; h < len; h++){
        texts_loc[h] = qatd_cpp_detect_hash_vector(texts[h], texts_loc[h], seq, id);
    }
    return(texts_loc);
}

/***R

toks_hash <- list(rep(1:10, 10), rep(5:15, 10))
loc <- qatd_cpp_structcopy_int_list(toks_hash)
loc <- qatd_cpp_detect_hash_list(toks_hash, loc, c(1, 2), 99)
loc <- qatd_cpp_detect_hash_list(toks_hash, loc, c(5, 6), 88)
loc

toks_hash <- rep(1:100, 1000)
microbenchmark::microbenchmark(
    qatd_cpp_detect_hash_vector(toks_hash, rep(0, 100000), c(3, 4), 9999),
    quanteda:::matchSequence(c(3, 4), toks_hash), times=100
)

*/
