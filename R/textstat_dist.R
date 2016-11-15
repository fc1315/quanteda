#' @title Distance matrix between documents and/or features 
#' 
#' @description These functions compute distance matrix between documents and/or features from a 
#' \code{\link{dfm}} and return a standard \code{\link[stats]{dist}} object.  
#'     
#' @slot selection character or character vector of document names or feature 
#'   labels from the dfm
#' @slot n the top \code{n} most similar items will be returned, sorted in 
#'   descending order.  If n is \code{NULL}, return all items.
#' @slot margin identifies the margin of the dfm on which similarity will be 
#'   computed:  \code{documents} for documents or \code{features} for word/term
#'   features.
#' @slot method a valid method for computing similarity from 
#'   \code{\link[proxy]{pr_DB}}
#' @slot sorted sort results in descending order if \code{TRUE}
#' @slot ordered whether a term appears before or after the target feature 
#'      are counted seperately
#' @slot normalize a deprecated argument retained (temporarily) for legacy 
#'   reasons.  If you want to compute similarity on a "normalized" dfm objects 
#'   (e.g. \code{x}), wrap it in \code{\link{weight}(x, "relFreq")}.
#' @slot digits decimal places to display similarity values
#' @slot tri whether the upper triangle of the symmetric \eqn{V \times V} matrix is recorded
#' @slot diag whether the diagonal of the distance matrix should be recorded
#' @seealso \link{dfm}
#' @export
#' @import methods
#' @docType class
#' @name dist-class
#' @keywords internal
setClass("dist",
         slots = c(selection = "character", n = "integer", margin = "character", 
                   method = "character", ordered = "logical", normalize = "logical", sorted = "logical",
                   digits = "integer", tri = "logical", diag = "logical"),
         prototype = list (selection = character(0), n = NULL,
                           margin = c("documents", "features"),
                           method = "correlation",
                           sorted = TRUE, normalize = FALSE, tri = FALSE, diag = FALSE)
         #contains = "Matrix"
         )

#' compute distance matrix between documents and/or features from a 
#' \code{\link{dfm}}. 

#' @examples
#' # create a dfm from inaugural addresses from Reagan onwards
#' presDfm <- dfm(subset(inaugCorpus, Year > 1980), ignoredFeatures = stopwords("english"),
#'                stem = TRUE)
#' 
#' # compute some document similarities
#' (tmp <- textstat_dist(presDfm, margin = "documents"))
#' # output as a matrix
#' as.matrix(tmp)
#' # for specific comparisons
#' textstat_dist(presDfm, "1985-Reagan", n = 5, margin = "documents")
#' textstat_dist(presDfm, c("2009-Obama" , "2013-Obama"), n = 5, margin = "documents")
#' textstat_dist(presDfm, c("2009-Obama" , "2013-Obama"), margin = "documents")
#' textstat_dist(presDfm, c("2009-Obama" , "2013-Obama"), margin = "documents", method = "cosine")
#' #textstat_dist(presDfm, "2005-Bush", margin = "documents", method = "eJaccard", sorted = FALSE)
#' 
#' # compute some term similarities
#' textstat_dist(presDfm, c("fair", "health", "terror"), method="cosine", margin = "features", 20)
#' 
#' @rdname dist-class
#' @export
setGeneric(name = "textstat_dist",
           signature = c("x", "selection", "n", "margin", "sorted", "normalize"),
           def = function(x, selection = character(0), n = NULL,
                          margin = c("documents", "features"),
                          method = "correlation",
                          sorted = TRUE, normalize = FALSE, tri = FALSE, diag = FALSE)
               {
                standardGeneric("textstat_dist")
               })


#' @rdname dist-class
#' @export
setMethod(f = "textstat_dist", 
          signature = signature("dfm", "ANY"),
          def = function(x, selection = character(0), n = NULL, 
                         margin = c("documents", "features"),
                         method = "correlation",
                         sorted = TRUE, normalize = FALSE, tri = TRUE, diag = FALSE ) {
              
              # value <- match.arg(value)
              
              if (normalize) {
                  warning("normalize is deprecated, not applied - use weight() instead")
                  # x <- weight(x, "relFreq")  # normalize by term freq.
              }
              
              margin <- match.arg(margin)
              if (margin == "features") {
                  items <- features(x)
              } else {
                  items <- docnames(x)
              }
              
              if (is.null(n) || n >= length(items))
                  #n <- length(items) - 1 # choose all features/docs if n is NULL
                  n <- length(items) 
              
              if (length(selection) != 0L) {
                  # retain only existing features or documents
                  selectIndex <- which(items %in% selection)
                  if (length(selectIndex)==0)
                      stop("no such documents or feature labels exist.")
                  
                  if (margin=="features") {
                      xSelect <- x[, selectIndex, drop=FALSE]
                  } else {
                      xSelect <- x[selectIndex, , drop=FALSE]
                  }
              } else xSelect <- NULL
              
              if (method == "cosine") {
                  result <- cosineSparse(x, xSelect, margin = ifelse(margin == "documents", 1, 2))
              } else if (method == "correlation") {
                  result <- correlationSparse(x, xSelect, margin = ifelse(margin == "documents", 1, 2))
              } else if (method == "Euclidean") {
                  result <- euclideanSparse(x, xSelect, margin = ifelse(margin == "documents", 1, 2))
              } else {
                  # use proxy::dist() for all other methods
                  result <- as.matrix(proxy::dist(as.matrix(x), as.matrix(xSelect), method = method, 
                                                        by_rows = ifelse(margin=="features", FALSE, TRUE)), diag = 1)
              }
              
              # rounding
              # if (!is.null(digits)) similmatrix <- round(similmatrix, digits)
              
              # convert NaNs to NA
              # similmatrix[is.nan(similmatrix)] <- NA
              
              # return the matrix if this was requested
              # return(similmatrix)
              
              # convert the matrix to a list of similarities
              # result <- lapply(seq_len(ncol(similmatrix)), function(i) similmatrix[, i])
              # names(result) <- if (!is.null(xSelect)) items[selectIndex] else if (margin == "documents") docnames(x) else features(x)
              
              # remove the element of each similarity vector equal to the item itself
              #tempseq <- seq_along(result)
              #names(tempseq) <- names(result)
              #result <- lapply( tempseq, function(i)
              #    result[[i]] <- result[[i]][-which(names(result[[i]]) == names(result)[i])] )
              
              # sort each element of the list and return only first n results if n not NULL
              if (sorted == TRUE)
                  # result <- lapply(result, sort, decreasing=TRUE, na.last = TRUE)
                  result[do.call(order, lapply(1:NCOL(result), function(i) result[, i])), ]

              # truncate to n if n is not NULL
              if (!is.null(n))
                  result <- head(result, n)

              # discard the upper diagonal if tri == TRUE
              if (tri)
                  result[upper.tri(result, diag = !diag)]<-0

              # create a new feature context matrix
              # result <- new("textstat_dist",selection = selection, n = n,
              #               margin = margin,
              #               method = method,
              #               sorted = sorted, normalize = normalize, tri = tri, diag = diag)
              class(result) <- "dist"
              # set the dimnames of result
              #dimnames(result) <- list(features = rownames(result), features = rownames(result))
              
              
              attr(result, "Size") <- nrow(result)
              if (!is.null(rownames(result))) attr(result, "Labels") <- rownames(result)
              attr(result, "Diag") <- diag
              attr(result, "Upper") <- !tri
              attr(result, "method") <- method
              attr(result, "call") <- match.call()
              attr(result, "dimnames") <- NULL
              result
          })


#' @rdname dist-class
#' @param ... unused
#' @export
as.list.dist <- function(x, ...) {
    
}



## code below based on assoc.R from the qlcMatrix package
## used Matrix::crossprod and Matrix::tcrossprod for sparse Matrix handling

# L2 norm
norm2 <- function(x,s) { drop(Matrix::crossprod(x^2, s)) ^ 0.5 }
# L1 norm
norm1 <- function(x,s) { drop(Matrix::crossprod(abs(x),s)) }

cosineSparse <- function(x, y = NULL, margin = 1, norm = norm2) {
    if (!(margin %in% 1:2)) stop("margin can only be 1 (rows) or 2 (columns)")
    if (margin == 1) x <- t(x)
    S <- rep(1, nrow(x))			
    N <- Matrix::Diagonal( x = norm(x, S)^-1 )
    x <- x %*% N
    if (!is.null(y)) {
        if (margin == 1) y <- t(y)
        N <- Matrix::Diagonal( x = match.fun(norm)(y, S)^-1 )
        y <- y %*% N
        return(as.matrix(Matrix::crossprod(x,y)))
    } else
        return(as.matrix(Matrix::crossprod(x)))
}

correlationSparse <- function(x, y = NULL, margin = 1) {
    if (!(margin %in% 1:2)) stop("margin can only be 1 (rows) or 2 (columns)")
    cpFun <- if (margin == 2) Matrix::crossprod else Matrix::tcrossprod
    tcpFun <- if (margin == 2) Matrix::tcrossprod else Matrix::crossprod
    marginSums <- if (margin == 2) colSums else rowSums
    
    n <- if (margin == 2) nrow(x) else ncol(x)
    muX <- if (margin == 2) colMeans(x) else rowMeans(x)
    
    if (!is.null(y)) {
        stopifnot(ifelse(margin == 2, nrow(x) == nrow(y), ncol(x) == ncol(y)))
        muY <- if (margin == 2) colMeans(y) else rowMeans(y)
        covmat <- (as.matrix(cpFun(x,y)) - n * tcrossprod(muX, muY)) / (n-1)
        sdvecX <- sqrt((marginSums(x^2) - n * muX^2) / (n-1))
        sdvecY <- sqrt((marginSums(y^2) - n * muY^2) / (n-1))
        cormat <- covmat / tcrossprod(sdvecX, sdvecY)
    } else {
        covmat <- ( as.matrix(cpFun(x)) - drop(n * tcrossprod(muX)) ) / (n-1)
        sdvec <- sqrt(diag(covmat))
        cormat <- covmat / tcrossprod(sdvec)
    }
    cormat
}

euclideanSparse <- function(x, y = NULL, margin = 1){
    if (!(margin %in% 1:2)) stop("margin can only be 1 (rows) or 2 (columns)")
    marginSums <- if (margin == 2) colSums else rowSums
    n <- if (margin == 2) ncol(x) else nrow(x)
    
    if (!is.null(y)) {
        stopifnot(ifelse(margin == 2, nrow(x) == nrow(y), ncol(x) == ncol(y)))
        an <- marginSums(x^2)
        bn <- marginSums(y^2)
        tmp <- matrix(rep(an, n), nrow = n) 
        tmp <-  tmp +  matrix(rep(bn, n), nrow = n, byrow=TRUE)
        eucmat <- sqrt( tmp - 2 * Matrix::tcrossprod(x, y) )
    } else {
        an <- marginSums(x^2)
        tmp <- matrix(rep(an, n), nrow = n) 
        tmp <-  tmp +  matrix(rep(an, n), nrow = n, byrow=TRUE)
        eucmat <- sqrt( tmp - 2 * as.matrix(Matrix::tcrossprod(x)))
    }
    eucmat
}
## JACCARD SPARSE:  See http://stackoverflow.com/questions/36220585/efficient-jaccard-similarity-documenttermmatrix

