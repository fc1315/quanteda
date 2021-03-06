# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

fcm_cpp <- function(texts, types, count, window, weights, ordered, tri, nvec) {
    .Call('quanteda_fcm_cpp', PACKAGE = 'quanteda', texts, types, count, window, weights, ordered, tri, nvec)
}

fcm_hash_cpp <- function(texts, n_types, count, window, weights, ordered, tri, nvec) {
    .Call('quanteda_fcm_hash_cpp', PACKAGE = 'quanteda', texts, n_types, count, window, weights, ordered, tri, nvec)
}

qatd_cpp_ngram_hashed_vector <- function(tokens, ns, skips) {
    .Call('quanteda_qatd_cpp_ngram_hashed_vector', PACKAGE = 'quanteda', tokens, ns, skips)
}

qatd_cpp_ngram_hashed_list <- function(texts, ns, skips) {
    .Call('quanteda_qatd_cpp_ngram_hashed_list', PACKAGE = 'quanteda', texts, ns, skips)
}

qatd_cpp_ngram_unhash_type <- function(ids_ngram, tokens, delim) {
    .Call('quanteda_qatd_cpp_ngram_unhash_type', PACKAGE = 'quanteda', ids_ngram, tokens, delim)
}

skipgramcpp <- function(tokens, ns, ks, delim) {
    .Call('quanteda_skipgramcpp', PACKAGE = 'quanteda', tokens, ns, ks, delim)
}

find_sequence_cppl <- function(texts, types, count_min, smooth, nested) {
    .Call('quanteda_find_sequence_cppl', PACKAGE = 'quanteda', texts, types, count_min, smooth, nested)
}

qatd_cpp_detect_hash_vector <- function(tokens_, tokens_loc_, seq_, id) {
    .Call('quanteda_qatd_cpp_detect_hash_vector', PACKAGE = 'quanteda', tokens_, tokens_loc_, seq_, id)
}

qatd_cpp_detect_hash_list <- function(texts_, texts_loc_, seq, id) {
    .Call('quanteda_qatd_cpp_detect_hash_list', PACKAGE = 'quanteda', texts_, texts_loc_, seq, id)
}

qatd_cpp_replace_hash_vector <- function(tokens_, seq_, id) {
    .Call('quanteda_qatd_cpp_replace_hash_vector', PACKAGE = 'quanteda', tokens_, seq_, id)
}

qatd_cpp_replace_hash_list <- function(texts_, flags, seq, id) {
    .Call('quanteda_qatd_cpp_replace_hash_list', PACKAGE = 'quanteda', texts_, flags, seq, id)
}

join_tokens_cpp <- function(tokens, tokens_join, delim) {
    invisible(.Call('quanteda_join_tokens_cpp', PACKAGE = 'quanteda', tokens, tokens_join, delim))
}

join_tokens_cppl <- function(texts, flags, tokens_join, delim) {
    invisible(.Call('quanteda_join_tokens_cppl', PACKAGE = 'quanteda', texts, flags, tokens_join, delim))
}

select_tokens_cppl <- function(texts, flags, types, remove, spacer) {
    invisible(.Call('quanteda_select_tokens_cppl', PACKAGE = 'quanteda', texts, flags, types, remove, spacer))
}

qatd_cpp_lookup_int_list <- function(texts_, texts_loc_, keys, id) {
    .Call('quanteda_qatd_cpp_lookup_int_list', PACKAGE = 'quanteda', texts_, texts_loc_, keys, id)
}

qatd_cpp_deepcopy <- function(x_) {
    .Call('quanteda_qatd_cpp_deepcopy', PACKAGE = 'quanteda', x_)
}

qatd_cpp_structcopy_int_list <- function(list_) {
    .Call('quanteda_qatd_cpp_structcopy_int_list', PACKAGE = 'quanteda', list_)
}

qatd_cpp_remove_chr_list <- function(list_, elem_remove) {
    .Call('quanteda_qatd_cpp_remove_chr_list', PACKAGE = 'quanteda', list_, elem_remove)
}

qatd_cpp_remove_int_list <- function(list_, elem_remove) {
    .Call('quanteda_qatd_cpp_remove_int_list', PACKAGE = 'quanteda', list_, elem_remove)
}

wordfishcpp <- function(wfm, dir, priors, tol, disp, dispfloor) {
    .Call('quanteda_wordfishcpp', PACKAGE = 'quanteda', wfm, dir, priors, tol, disp, dispfloor)
}

