#' Create ngrams and skipgrams
#' 
#' Create a set of ngrams (tokens in sequence) from character vectors or
#' tokenized text objects, with an optional skip argument to form skipgrams. 
#' Both the ngram length and the skip lengths take vectors of arguments to form 
#' multiple lengths or skips in one pass.  \code{ngrams()} is implemented in C++
#' for efficiency.
#' @author Kohei Watanabe and Ken Benoit
#' @return a tokenizedTexts object consisting a list of character vectors of 
#'   ngrams, one list element per text, or a character vector if called on a 
#'   simple character vector
#' @param x a tokenizedText object or a character vector of tokens
#' @param n integer vector specifying the number of elements to be concatenated 
#'   in each ngram
#' @param skip integer vector specifying the adjacency skip size for tokens 
#'   forming the ngrams, default is 0 for only immediately neighbouring words. 
#'   For \code{skipgrams}, \code{skip} can be a vector of integers, as the 
#'   "classic" approach to forming skip-grams is to set skip = \eqn{k} where
#'   \eqn{k} is the distance for which \eqn{k} or fewer skips are used to construct
#'   the \eqn{n}-gram.  Thus a "4-skip-n-gram" defined as \code{skip = 0:4}
#'   produces results that include 4 skips, 3 skips, 2 skips, 1 skip, and 0
#'   skips (where 0 skips are typical n-grams formed from adjacent words).  See
#'   Guthrie et al (2006).
#' @param concatenator character for combining words, default is \code{_} 
#'   (underscore) character
#' @param ... not used
#' @export
#' @examples
#' # ngrams
#' tokens <- tokens("the quick brown fox jumped over the lazy dog.", 
#'                    removePunct = TRUE, simplify = TRUE)
#' ngrams(tokens, n = 1:3)
#' ngrams(tokens, n = c(2,4), concatenator = " ")
#' ngrams(tokens, n = c(2,4), skip = 1, concatenator = " ")
#' 
#' # skipgrams
ngrams <- function(x, ...) {
    UseMethod("ngrams")
}

#' @rdname ngrams
#' @importFrom stats complete.cases
#' @export
ngrams.character <- function(x, n = 2L, skip = 0L, concatenator = "_", ...) {
    # trap condition where a "text" is a single NA
    if (is.na(x[1]) && length(x)==1) return(NULL)
    if (any(stringi::stri_detect_fixed(x, " ")) & concatenator != " ")
        warning("whitespace detected: you may need to run tokens() first")
    if (length(x) < min(n)) return(NULL)
    if (identical(as.integer(n), 1L)) {
        if (!identical(as.integer(skip), 0L))
            warning("skip argument ignored for n = 1")
        return(x)
    }
    skipgramcpp(x, n, skip + 1, concatenator)
}

#' @rdname ngrams
#' @examples 
#' # with "classic" tokenized objects
#' txt <- c("a b c d e", "c d e f g")
#' toks <- tokenize(txt)
#' ngrams(toks, n = 2:3)
#' @export
ngrams.tokenizedTexts <- function(x, ...) {
    as.tokenizedTexts(ngrams(as.tokens(x), ...))
}

#' @rdname ngrams
#' @examples 
#' txt <- c("a b c d e", "c d e f g")
#' toks_hashed <- tokens(txt)
#' ngrams(toks_hashed, n = 2:3)
#' @export
ngrams.tokens <- function(x, n = 2L, skip = 0L, concatenator = "_", ...) {
    attrs_orig <- attributes(x)
    if (any(n <= 0)) stop("ngram length has to be greater than zero")
    # Generate ngrams
    res <- qatd_cpp_ngram_hashed_list(x, n, skip + 1)
    
    # Make character tokens of ngrams
    ngram_types <- qatd_cpp_ngram_unhash_type(res$id_unigram, types(x), concatenator)
    ngramsResult <- res$text
    names(ngramsResult) <- names(x)
    attributes(ngramsResult) <- attrs_orig
    types(ngramsResult) <- ngram_types
    attr(ngramsResult, "ngrams") <- as.integer(n)
    attr(ngramsResult, "skip") <- as.integer(skip)
    attr(ngramsResult, "concatenator") <- concatenator
    ngramsResult
}

#' @rdname ngrams
#' @details Normally, \code{\link{ngrams}} will be called through 
#'   \code{\link{tokenize}}, but these functions are also exported in case a 
#'   user wants to perform lower-level ngram construction on tokenized texts.
#'   
#'   \code{\link{skipgrams}} is a wrapper to \code{\link{ngrams}} that requires 
#'   arguments to be supplied for both \code{n} and \code{skip}.  For \eqn{k}-skip 
#'   skipgrams, set \code{skip} to \code{0:}\eqn{k}, in order to conform to the
#'   definition of skip-grams found in Guthrie et al (2006): A \eqn{k} skip-gram is
#'   an ngram which is a superset of all ngrams and each \eqn{(k-i)} skipgram until
#'   \eqn{(k-i)==0} (which includes 0 skip-grams).
#' @export
#' @references 
#' \href{http://homepages.inf.ed.ac.uk/ballison/pdf/lrec_skipgrams.pdf}{Guthrie,
#' D., B. Allison, W. Liu, and L. Guthrie. 2006. "A Closer Look at Skip-Gram 
#' Modelling."}
#' @importFrom utils combn
#' @examples 
#' tokens <- tokens("insurgents killed in ongoing fighting", simplify = TRUE)
#' skipgrams(tokens, n = 2, skip = 0:1, concatenator = " ") 
#' skipgrams(tokens, n = 2, skip = 0:2, concatenator = " ") 
#' skipgrams(tokens, n = 3, skip = 0:2, concatenator = " ")   
skipgrams <- function(x, ...) UseMethod("skipgrams")

#' @rdname ngrams
#' @export
skipgrams.character <- function(x, n, skip, concatenator="_", ...)
    ngrams(x, n, skip, concatenator)

#' @rdname ngrams
#' @export
skipgrams.tokenizedTexts <- function(x, n, skip, concatenator="_", ...)
    ngrams(x, n, skip, concatenator, ...)

#' @rdname ngrams
#' @export
skipgrams.tokens <- function(x, n, skip, concatenator="_", ...)
    ngrams(x, n, skip, concatenator, ...)
