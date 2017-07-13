#' Fan and Vercauteren encryption scheme
#'
#' The Fan and Vercauteren scheme is implemented in this package.
#'
#' Description of the scheme.
#'
#' @name FandV
#'
#' @usage
#' p <- pars("FandV")
#'
#' @examples
#' # Benchmark the performance of the scheme
#' library(microbenchmark)
#' p <- pars("FandV")
#' microbenchmark({ keys <- keygen(p) }, unit="ms")
#' microbenchmark({ ct1 <- enc(keys$pk, 2) }, unit="ms")
#' ct2 <- enc(keys$pk, 3)
#' microbenchmark({ ct1 + ct2 }, unit="ms")
#' microbenchmark({ ct1 * ct2 }, unit="ms")
#' microbenchmark({ dec(keys$sk, ct1) }, unit="ms")
#'
NULL

loadModule("FandV", TRUE)

rlkLocker <- 42 # Setup var in the namespace of the package ...
# http://stackoverflow.com/questions/18151619/operator-overloading-in-r-reference-classes
# evalqOnLoad used in package RcppBDT
evalqOnLoad({
  rlkLocker <<- new(FandV_rlk_locker) # ... then overwrite in evalqOnLoad once Rcpp modules loaded

  ##### Missing S4 generics #####
  setGeneric("diag")
  setGeneric("diag<-")
  setGeneric("t")
  setGeneric("dimnames<-")
  setGeneric("rowSums")
  setGeneric("colSums")
  setGeneric("crossprod")
  setGeneric("tcrossprod")

  ##### Single ciphertexts #####
  setMethod("+", c("Rcpp_FandV_ct", "Rcpp_FandV_ct"), function(e1, e2) {
    ct <- e1$add(e2)
    # Prepare return result
    attr(ct, "FHEt") <- "ct"
    attr(ct, "FHEs") <- "FandV"
    ct
  })
  setMethod("-", c("Rcpp_FandV_ct", "Rcpp_FandV_ct"), function(e1, e2) {
    ct <- e1$sub(e2)
    # Prepare return result
    attr(ct, "FHEt") <- "ct"
    attr(ct, "FHEs") <- "FandV"
    ct
  })
  setMethod("*", c("Rcpp_FandV_ct", "Rcpp_FandV_ct"), function(e1, e2) {
    ct <- e1$mul(e2)
    # Prepare return result
    attr(ct, "FHEt") <- "ct"
    attr(ct, "FHEs") <- "FandV"
    ct
  })
  setMethod("rep", signature(x="Rcpp_FandV_ct"), function(x, ...) {
    idx <- rep(1, ...)

    res <- new(FandV_ct_vec)
    for(i in idx) {
      res$push(x)
    }

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })

  ##### Vectors of ciphertexts #####
  ### TODO: diff
  setMethod("c", signature(x="Rcpp_FandV_ct"), function (x, ..., recursive = FALSE) {
    res <- new(FandV_ct_vec)
    res$push(x)

    if(!missing(...)) {
      args <- list(...)
      for(i in 1:length(args)) {
        if(class(args[[i]])=="Rcpp_FandV_ct") {
          res$push(args[[i]])
        } else if(class(args[[i]])=="Rcpp_FandV_ct_vec") {
          res$pushvec(args[[i]])
        } else {
          stop("only Fan and Vercauteren ciphertexts or ciphertext vectors can be concatenated")
        }
      }
    }

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("c", "Rcpp_FandV_ct_vec", function (x, ..., recursive = FALSE) {
    res <- new(FandV_ct_vec)
    res$pushvec(x)

    args <- list(...)
    for(i in 1:length(args)) {
      if(class(args[[i]])=="Rcpp_FandV_ct") {
        res$push(args[[i]])
      } else if(class(args[[i]])=="Rcpp_FandV_ct_vec") {
        res$pushvec(args[[i]])
      } else {
        stop("only Fan and Vercauteren ciphertexts or ciphertext vectors can be concatenated")
      }
    }

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("[", "Rcpp_FandV_ct_vec", function(x, i, j, ..., drop=TRUE) {
    i <- as.integer(i)
    if(max(abs(i))>x$size()) {
      stop("out of bounds")
    }
    if(min(i) < 0 && max(i) > 0) {
      stop("only 0's may be mixed with negative subscripts")
    }
    if(min(i) < 0)
      res <- x$without(-sort(i)-1)
    else
      res <- x$subset(i-1)
    if(res$size() == 1) {
      ct <- res$get(0)
      attr(ct, "FHEt") <- "ct"
      attr(ct, "FHEs") <- "FandV"
      return(ct)
    }
    if(res$size() == 0) {
      return(NULL)
    }

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_vec", value="Rcpp_FandV_ct"), function (x, i, j, ..., value) {
    i <- as.integer(i)
    idx <- (1:x$size())[i]
    if(any(is.na(idx))) {
      stop("out of bounds")
    }
    for(j in idx) {
      x$set(idx-1, value)
    }

    attr(x, "FHEt") <- "ctvec"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_vec", value="Rcpp_FandV_ct_vec"), function (x, i, j, ..., value) {
    i <- as.integer(i)
    idx <- (1:x$size())[i]
    if(any(is.na(idx))) {
      stop("out of bounds")
    }
    if(length(i)!=length(value)) {
      warning("number of items to replace is not a multiple of replacement length.")
    }
    for(j in 1:length(idx)) {
      x$set(idx[j]-1, value$get((j-1)%%value$size()))
    }

    attr(x, "FHEt") <- "ctvec"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_vec"), function (x, i, j, ..., value) {
    stop("only a ciphertext can be assigned to this vector")
  })
  setMethod("+", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    if(e1$size()%%e2$size()!=0 && e2$size()%%e1$size()!=0) {
      stop("longer object length is not a multiple of shorter object length")
    }
    res <- e1$add(e2)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("-", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    if(e1$size()%%e2$size()!=0 && e2$size()%%e1$size()!=0) {
      stop("longer object length is not a multiple of shorter object length")
    }
    res <- e1$sub(e2)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    if(e1$size()%%e2$size()!=0 && e2$size()%%e1$size()!=0) {
      warning("longer object length is not a multiple of shorter object length")
    }
    res <- e1$mulParallel(e2)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  # This is weird.  %*% doesn't support S4 method dispatch.  I think this is because
  # this pkg 'Depend's on gmp and for some reason they force S3 dispatch on %*%
  # See gmp package source: gmp/R/matrix-prods.R, line 47 (top is if(FALSE)'ed out)
#   setMethod("%*%", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct_vec"), function(x, y) {
#     if(x$size()!=y$size()) {
#       stop("non-conformable arguments")
#     }
#     res <- x$innerprod(y)
#
#     attr(res, "FHEt") <- "ct"
#     attr(res, "FHEs") <- "FandV"
#     res
#   }
  setMethod("+", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct"), function(e1, e2) {
    res <- e1$addct(e2)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("+", c("Rcpp_FandV_ct", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    res <- e2$addct(e1)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("-", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct"), function(e1, e2) {
    res <- e1$subct(e2, FALSE)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("-", c("Rcpp_FandV_ct", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    res <- e2$subct(e1, TRUE)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct"), function(e1, e2) {
    res <- e1$mulctParallel(e2)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    res <- e2$mulctParallel(e1)

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("rep", signature(x="Rcpp_FandV_ct_vec"), function(x, ...) {
    idx <- rep(1:length(x), ...)

    res <- new(FandV_ct_vec)
    for(i in idx) {
      res$push(x[i])
    }

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("sum", c("Rcpp_FandV_ct_vec", "logical"), function(x, na.rm) {
    if(x$size() < 40) {
      res <- x$sumSerial()
    } else if(defaultNumThreads()*20>x$size()) {
      setThreadOptions(x$size()%/%20)
      res <- x$sumParallel()
      setThreadOptions()
    } else {
      res <- x$sumParallel()
    }

    attr(res, "FHEt") <- "ct"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("prod", c("Rcpp_FandV_ct_vec", "logical"), function(x, na.rm) {
    res <- x$prodParallel()

    attr(res, "FHEt") <- "ct"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("length", signature(x="Rcpp_FandV_ct_vec"), function(x) {
    return(x$size())
  })
  setMethod("diag", signature(x="Rcpp_FandV_ct_vec", nrow="missing"), function(x, ncol) { diag(x, length(x), ncol) })
  setMethod("diag", signature(x="Rcpp_FandV_ct_vec"), function(x, nrow, ncol) {
    if(missing(nrow) && missing(ncol))
      nrow <- ncol <- length(x)
    if(missing(ncol) && !missing(nrow))
      ncol <- nrow

    tmp <- new(FandV_ct, x[1]$p, rlkLocker, x[1]$rlki)

    res <- matrix(tmp, nrow, ncol)
    for(i in 0:(min(c(nrow,ncol))-1)) {
      res$set(i, i, x[(i%%length(x))+1])
    }

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })

  ##### Matrices of ciphertexts #####
  # TODO: norm (those that can be done), diff
  # gmp package again overrides matrix and makes it S3 dispatch
#   setMethod("matrix", "Rcpp_FandV_ct_vec", function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, ...) {
#   })
  # TODO: probably not that useful, but base R will allow a matrix to concatenate
  #       with a scalar vector and just coerces the matrix to a vector columnwise
#   setMethod("c", "Rcpp_FandV_ct_mat", function (x, ..., recursive = FALSE) {
#   })
  setMethod("[", signature(x="Rcpp_FandV_ct_mat", i="missing", j="missing", drop="ANY"), function(x, i, j, ..., drop=TRUE) x)
  setMethod("[", signature(x="Rcpp_FandV_ct_mat", i="numeric", j="missing", drop="ANY"), function(x, i, j, ..., drop=TRUE) {
    i <- as.integer(i)
    if(max(abs(i))>x$size()) {
      stop("out of bounds")
    }
    if(min(i) < 0 && max(i) > 0) {
      stop("only 0's may be mixed with negative subscripts")
    }
    if(nargs() == 2) { # eg ct[1] or ct[1:4] etc
      tmp <- matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[i]
      res <- x$subsetV(as.vector(tmp))
      attr(res, "FHEt") <- "ctvec"
      attr(res, "FHEs") <- "FandV"
      if(res$size() == 1) {
        ct <- res$get(0)
        attr(ct, "FHEt") <- "ct"
        attr(ct, "FHEs") <- "FandV"
        return(ct)
      }
      if(res$size() == 0) {
        return(NULL)
      }

      return(res)
    } else { # eg ct[1,] or ct[1:4,]
      tmp <- matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[i,,drop=drop]
      if(is.matrix(tmp)) {
        res <- x$subset(as.vector(tmp), nrow(tmp), ncol(tmp))
        attr(res, "FHEt") <- "ctmat"
        attr(res, "FHEs") <- "FandV"
      } else {
        res <- x$subsetV(as.vector(tmp))
        attr(res, "FHEt") <- "ctvec"
        attr(res, "FHEs") <- "FandV"
      }
      if(res$size() == 1) {
        ct <- res$get(0)
        attr(ct, "FHEt") <- "ct"
        attr(ct, "FHEs") <- "FandV"
        return(ct)
      }
      if(res$size() == 0) {
        return(NULL)
      }

      return(res)
    }
  })
  setMethod("[", signature(x="Rcpp_FandV_ct_mat", i="missing", j="numeric", drop="ANY"), function(x, i, j, ..., drop=TRUE) {
    # eg ct[,1] or ct[,1:4]
    j <- as.integer(j)
    if(max(abs(j))>x$ncol) {
      stop("out of bounds")
    }
    if(min(j) < 0 && max(j) > 0) {
      stop("only 0's may be mixed with negative subscripts")
    }
    tmp <- matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[,j,drop=drop]
    if(is.matrix(tmp)) {
      res <- x$subset(as.vector(tmp), nrow(tmp), ncol(tmp))
      attr(res, "FHEt") <- "ctmat"
      attr(res, "FHEs") <- "FandV"
    } else {
      res <- x$subsetV(as.vector(tmp))
      attr(res, "FHEt") <- "ctvec"
      attr(res, "FHEs") <- "FandV"
    }
    if(res$size() == 1) {
      ct <- res$get(0)
      attr(ct, "FHEt") <- "ct"
      attr(ct, "FHEs") <- "FandV"
      return(ct)
    }
    if(res$size() == 0) {
      return(NULL)
    }

    return(res)
  })
  setMethod("[", signature(x="Rcpp_FandV_ct_mat", i="numeric", j="numeric", drop="ANY"), function(x, i, j, ..., drop=TRUE) {
    i <- as.integer(i)
    j <- as.integer(j)
    if(max(abs(i))>x$nrow || max(abs(j))>x$ncol) {
      stop("out of bounds")
    }
    if((min(i) < 0 && max(i) > 0) || (min(j) < 0 && max(j) > 0)) {
      stop("only 0's may be mixed with negative subscripts")
    }
    tmp <- matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[i,j,drop=drop]
    if(is.matrix(tmp)) {
      res <- x$subset(as.vector(tmp), nrow(tmp), ncol(tmp))
      attr(res, "FHEt") <- "ctmat"
      attr(res, "FHEs") <- "FandV"
    } else {
      res <- x$subsetV(as.vector(tmp))
      attr(res, "FHEt") <- "ctvec"
      attr(res, "FHEs") <- "FandV"
    }
    if(res$size() == 1) {
      ct <- res$get(0)
      attr(ct, "FHEt") <- "ct"
      attr(ct, "FHEs") <- "FandV"
      return(ct)
    }
    if(res$size() == 0) {
      return(NULL)
    }

    return(res)
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat", i="numeric", j="numeric", value="Rcpp_FandV_ct"), function (x, i, j, ..., value) {
    i <- as.integer(i)
    j <- as.integer(j)
    if(length(i) > 1 || length(j) > 1)
      stop("only single element assignment currently supported for FandV ciphertext vectors")
    if((min(i) < 1 && max(i) > x$nrow) || (min(j) < 1 && max(j) > x$ncol)) {
      stop("out of bounds")
    }
    x$set(i-1, j-1, value)

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat", i="missing", j="numeric", value="Rcpp_FandV_ct_vec"), function (x, i, j, ..., value) {
    j <- as.integer(j)
    tmp <- c(matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[,j])
    if(any(is.na(tmp))) {
      stop("out of bounds")
    }
    if(length(tmp)%%length(value)!=0) {
      stop("number of items to replace is not a multiple of replacement length.")
    }
    for(k in 1:length(tmp)) {
      x$setelt(tmp[k], value[(k-1)%%length(value)+1])
    }

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat", i="numeric", j="missing", value="Rcpp_FandV_ct_vec"), function (x, i, j, ..., value) {
    i <- as.integer(i)
    tmp <- c(matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[i,])
    if(any(is.na(tmp))) {
      stop("out of bounds")
    }
    if(length(tmp)%%length(value)!=0) {
      stop("number of items to replace is not a multiple of replacement length.")
    }
    for(k in 1:length(tmp)) {
      x$setelt(tmp[k], value[(k-1)%%length(value)+1])
    }

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat", i="numeric", j="numeric", value="Rcpp_FandV_ct_vec"), function (x, i, j, ..., value) {
    i <- as.integer(i)
    j <- as.integer(j)
    tmp <- c(matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)[i,j])
    if(any(is.na(tmp))) {
      stop("out of bounds")
    }
    if(length(tmp)%%length(value)!=0) {
      stop("number of items to replace is not a multiple of replacement length.")
    }
    for(k in 1:length(tmp)) {
      x$setelt(tmp[k], value[(k-1)%%length(value)+1])
    }

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat", value="Rcpp_FandV_ct"), function (x, i, j, ..., value) {
    stop("only single element assignment currently supported for FandV ciphertext vectors")
  })
  setMethod("[<-", signature(x="Rcpp_FandV_ct_mat"), function (x, i, j, ..., value) {
    stop("only a ciphertext can be assigned to this vector")
  })
  setMethod("+", signature(e1="Rcpp_FandV_ct_mat", e2="Rcpp_FandV_ct_mat"), function(e1, e2) {
    if(e1$nrow!=e2$nrow || e2$ncol!=e2$ncol) {
      stop("non-conformable matrix sizes")
    }
    res <- e1$add(e2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", signature(e1="Rcpp_FandV_ct_mat", e2="Rcpp_FandV_ct_mat"), function(e1, e2) {
    if(e1$nrow!=e2$nrow || e2$ncol!=e2$ncol) {
      stop("non-conformable matrix sizes")
    }
    res <- e1$mul(e2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  # This is weird.  %*% doesn't support S4 method dispatch.  I think this is because
  # this pkg 'Depend's on gmp and for some reason they force S3 dispatch on %*%
  # See gmp package source: gmp/R/matrix-prods.R, line 47 (top is if(FALSE)'ed out)
#   setMethod("%*%", c("Rcpp_FandV_ct_mat", "Rcpp_FandV_ct_mat"), function(x, y) {})
  setMethod("+", c("Rcpp_FandV_ct_mat", "Rcpp_FandV_ct"), function(e1, e2) {
    res <- e1$addct(e2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("+", c("Rcpp_FandV_ct", "Rcpp_FandV_ct_mat"), function(e1, e2) {
    res <- e2$addct(e1)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct_mat", "Rcpp_FandV_ct"), function(e1, e2) {
    res <- e1$mulctParallel(e2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct", "Rcpp_FandV_ct_mat"), function(e1, e2) {
    res <- e2$mulctParallel(e1)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct_mat", "Rcpp_FandV_ct_vec"), function(e1, e2) {
    if(length(e2) > nrow(e1)*ncol(e1))
      stop("dims [product ", nrow(e1)*ncol(e1),"] do not match the length of object [", length(e2), "]")
    if((ncol(e1)*nrow(e1))%%length(e2)!=0)
      warning("longer object length is not a multiple of shorter object length")
    res <- e1$mulctvecParallel(e2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("*", c("Rcpp_FandV_ct_vec", "Rcpp_FandV_ct_mat"), function(e1, e2) {
    if(length(e1) > nrow(e2)*ncol(e2))
      stop("dims [product ", nrow(e2)*ncol(e2),"] do not match the length of object [", length(e1), "]")
    if((ncol(e2)*nrow(e2))%%length(e1)!=0)
      warning("longer object length is not a multiple of shorter object length")
    res <- e2$mulctvecParallel(e1)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_mat"), function(x, y) {
    res <- x$TmatmulParallel(y)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(nrow(x)!=length(y))
      stop("non-conformable arguments")
    y2 <- cbind(y)
    res <- x$TmatmulParallel(y2)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(nrow(y)!=length(x))
      stop("non-conformable arguments")
    x2 <- cbind(x)
    res <- x2$TmatmulParallel(y)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(length(x)!=length(y))
      stop("non-conformable arguments")
    crossprod(cbind(x), cbind(y))
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_mat", y="missing"), function(x, y) {
    crossprod(x, x)
  })
  setMethod("crossprod", signature(x="Rcpp_FandV_ct_vec", y="missing"), function(x, y) {
    crossprod(x, x)
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_mat"), function(x, y) {
    res <- x$matmulTParallel(y)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_vec"), function(x, y) {
    stop("non-conformable arguments") # It seems R will always error out no matter the vector size of the second argument in tcrossprod???  Assumes all vectors are column here, but elsewhere doesn't?
#     if(ncol(x)!=length(y))
#       stop("non-conformable arguments")
#     y2 <- rbind(y)
#     res <- x$matmulTParallel(y2)
#
#     attr(res, "FHEt") <- "ctmat"
#     attr(res, "FHEs") <- "FandV"
#     res
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(ncol(y)!=length(x))
      stop("non-conformable arguments")
    x2 <- rbind(x)
    res <- x2$matmulTParallel(y)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_vec"), function(x, y) {
    tcrossprod(cbind(x), cbind(y))
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_mat", y="missing"), function(x, y) {
    tcrossprod(x, x)
  })
  setMethod("tcrossprod", signature(x="Rcpp_FandV_ct_vec", y="missing"), function(x, y) {
    tcrossprod(x, x)
  })
  setMethod("dim", signature(x="Rcpp_FandV_ct_mat"), function(x) {
    c(x$nrow, x$ncol)
  })
  setMethod("length", signature(x="Rcpp_FandV_ct_mat"), function(x) {
    x$nrow*x$ncol
  })
  setMethod("diag", signature(x="Rcpp_FandV_ct_mat"), function(x, nrow, ncol) {
    if(!missing(nrow) || !missing(ncol))
      stop("'nrow' or 'ncol' cannot be specified when 'x' is a cipher text matrix")

    tmp <- matrix(0:(x$size()-1), nrow=x$nrow, ncol=x$ncol)
    res <- x$subsetV(as.vector(diag(tmp)))

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("diag<-", signature(x="Rcpp_FandV_ct_mat", value="Rcpp_FandV_ct_vec"), function(x, value) {
    if(length(value) != min(x$nrow, x$ncol))
      stop("replacement diagonal has wrong length")

    for(i in 0:(min(x$nrow, x$ncol)-1)) {
      x$set(i, i, value[i+1])
    }

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("diag<-", signature(x="Rcpp_FandV_ct_mat", value="Rcpp_FandV_ct"), function(x, value) {
    for(i in 0:(min(x$nrow, x$ncol)-1)) {
      x$set(i, i, value)
    }

    attr(x, "FHEt") <- "ctmat"
    attr(x, "FHEs") <- "FandV"
    x
  })
  setMethod("t", signature(x="Rcpp_FandV_ct_mat"), function(x) {
    res <- x$t()

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("rowSums", signature(x="Rcpp_FandV_ct_mat"), function(x, ...) {
    res <- x$rowSumsParallel()

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("colSums", signature(x="Rcpp_FandV_ct_mat"), function(x, ...) {
    res <- x$colSumsParallel()

    attr(res, "FHEt") <- "ctvec"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("dimnames<-", signature(x="Rcpp_FandV_ct_mat", value="ANY"), function(x, value) { x }) # Dummy so that rbind/cbind works
  ### RBIND
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_mat", y="missing"), function(x, y) { x })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_mat", y="NULL"), function(x, y) { x })
  setMethod("rbind2", signature(x="NULL", y="Rcpp_FandV_ct_mat"), function(x, y) { y })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_vec", y="missing"), function(x, y) {
    res <- new(FandV_ct_mat)
    res$setmatrix(x, 1, length(x), TRUE)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_vec", y="NULL"), function(x, y) { rbind2(x) })
  setMethod("rbind2", signature(x="NULL", y="Rcpp_FandV_ct_vec"), function(x, y) { rbind2(y) })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct", y="missing"), function(x, y) {
    res <- new(FandV_ct_mat)
    res$reset(x, 1, 1)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct", y="NULL"), function(x, y) { rbind2(x) })
  setMethod("rbind2", signature(x="NULL", y="Rcpp_FandV_ct"), function(x, y) { rbind2(y) })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(ncol(x)!=ncol(y))
      stop("number of columns of matrices must match")

    nrx <- nrow(x)
    nry <- nrow(y)
    nc <- ncol(x)
    # Create new matrix of the right size and fill with holding data
    res <- new(FandV_ct_mat)
    res$reset(x[1,1], nrx+nry, nc)

    ## FILL
    for(i in 0:(nrx-1)) {
      for(j in 0:(nc-1)) {
        res$set(i, j, x$get(i + j*nrx))
      }
    }
    for(i in nrx:(nrx+nry-1)) {
      for(j in 0:(nc-1)) {
        res$set(i, j, y$get(i-nrx + j*nry))
      }
    }

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(ncol(x)%%length(y)!=0)
      warning("number of columns of result is not a multiple of vector length")

    # Make y into a matrix and then we can be lazy and use the above method
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, 1, ncol(x), TRUE)

    rbind2(x, y2)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(ncol(y)%%length(x)!=0)
      warning("number of columns of result is not a multiple of vector length")

    # Make y into a matrix and then we can be lazy and use the above method
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, 1, ncol(y), TRUE)

    rbind2(x2, y)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct"), function(x, y) {
    # Make y into a matrix and then we can be lazy and use the above method
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), 1, ncol(x), TRUE)

    rbind2(x, y2)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct_mat"), function(x, y) {
    # Make y into a matrix and then we can be lazy and use the above method
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), 1, ncol(y), TRUE)

    rbind2(x2, y)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, 1, length(x), TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), 1, length(x), TRUE)
    rbind2(x2, y2)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct_vec"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), 1, length(y), TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, 1, length(y), TRUE)
    rbind2(x2, y2)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(length(x) > length(y))
      nc <- length(x)
    else
      nc <- length(y)
    if(nc%%length(x)!=0 || nc%%length(y)!=0)
      warning("number of columns of result is not a multiple of vector length")

    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, 1, nc, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, 1, nc, TRUE)
    rbind2(x2, y2)
  })
  setMethod("rbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), 1, 1, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), 1, 1, TRUE)
    rbind2(x2, y2)
  })
  ### CBIND
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_mat", y="missing"), function(x, y) { x })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_mat", y="NULL"), function(x, y) { x })
  setMethod("cbind2", signature(x="NULL", y="Rcpp_FandV_ct_mat"), function(x, y) { y })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_vec", y="missing"), function(x, y) {
    res <- new(FandV_ct_mat)
    res$setmatrix(x, length(x), 1, TRUE)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_vec", y="NULL"), function(x, y) { cbind2(x) })
  setMethod("cbind2", signature(x="NULL", y="Rcpp_FandV_ct_vec"), function(x, y) { cbind2(y) })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct", y="missing"), function(x, y) {
    res <- new(FandV_ct_mat)
    res$reset(x, 1, 1)

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct", y="NULL"), function(x, y) { cbind2(x) })
  setMethod("cbind2", signature(x="NULL", y="Rcpp_FandV_ct"), function(x, y) { cbind2(y) })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(nrow(x)!=nrow(y))
      stop("number of rows of matrices must match")

    ncx <- ncol(x)
    ncy <- ncol(y)
    nr <- nrow(x)
    # Create new matrix of the right size and fill with holding data
    res <- new(FandV_ct_mat)
    res$reset(x[1,1], nr, ncx+ncy)

    ## FILL
    for(i in 0:(nr-1)) {
      for(j in 0:(ncx-1)) {
        res$set(i, j, x$get(i + j*nr))
      }
    }
    for(i in 0:(nr-1)) {
      for(j in ncx:(ncx+ncy-1)) {
        res$set(i, j, y$get(i + (j-ncx)*nr))
      }
    }

    attr(res, "FHEt") <- "ctmat"
    attr(res, "FHEs") <- "FandV"
    res
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(nrow(x)%%length(y)!=0)
      warning("number of rows of result is not a multiple of vector length")

    # Make y into a matrix and then we can be lazy and use the above method
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, nrow(x), 1, TRUE)

    cbind2(x, y2)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_mat"), function(x, y) {
    if(nrow(y)%%length(x)!=0)
      warning("number of rows of result is not a multiple of vector length")

    # Make y into a matrix and then we can be lazy and use the above method
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, nrow(y), 1, TRUE)

    cbind2(x2, y)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_mat", y="Rcpp_FandV_ct"), function(x, y) {
    # Make y into a matrix and then we can be lazy and use the above method
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), nrow(x), 1, TRUE)

    cbind2(x, y2)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct_mat"), function(x, y) {
    # Make y into a matrix and then we can be lazy and use the above method
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), nrow(y), 1, TRUE)

    cbind2(x2, y)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, length(x), 1, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), length(x), 1, TRUE)
    cbind2(x2, y2)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct_vec"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), length(y), 1, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, length(y), 1, TRUE)
    cbind2(x2, y2)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct_vec", y="Rcpp_FandV_ct_vec"), function(x, y) {
    if(length(x) > length(y))
      nr <- length(x)
    else
      nr <- length(y)
    if(nr%%length(x)!=0 || nr%%length(y)!=0)
      warning("number of rows of result is not a multiple of vector length")

    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(x, nr, 1, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(y, nr, 1, TRUE)
    cbind2(x2, y2)
  })
  setMethod("cbind2", signature(x="Rcpp_FandV_ct", y="Rcpp_FandV_ct"), function(x, y) {
    # Make both matrices
    x2 <- new(FandV_ct_mat)
    x2$setmatrix(c(x), 1, 1, TRUE)
    y2 <- new(FandV_ct_mat)
    y2$setmatrix(c(y), 1, 1, TRUE)
    cbind2(x2, y2)
  })
  
})

matrix.Rcpp_FandV_ct_vec <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE, ...) {
  # Similar logic to r-source/src/main/array.c, do_matrix function from R
  if(missing(nrow) && missing(ncol)) {
    nrow <- length(data)
  } else if(missing(nrow)) {
    nrow <- ceiling(length(data)/ncol)
  } else if(missing(ncol)) {
    ncol <- ceiling(length(data)/nrow)
  }
  if((nrow*ncol)%%length(data)!=0) {
    stop("data length [", length(data), "] is not a sub-multiple or multiple of the number of rows/cols")
  }

  res <- new(FandV_ct_mat)
  res$setmatrix(data, nrow, ncol, byrow)

  attr(res, "FHEt") <- "ctmat"
  attr(res, "FHEs") <- "FandV"
  res
}
matrix.Rcpp_FandV_ct <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE, ...) {
  x <- new(FandV_ct_vec)
  x$push(data)
  matrix(x, nrow, ncol, byrow, ...)
}

# See above for why this is here
`%*%.Rcpp_FandV_ct_vec` <- function(x, y) {
  if(class(y) == "Rcpp_FandV_ct_mat")
    return(crossprod(x, y))
  if(class(y) != "Rcpp_FandV_ct_vec")
    stop("requires cipher text matrix/vector arguments")
  if(x$size()!=y$size())
    stop("non-conformable arguments")
  res <- x$innerprod(y)

  attr(res, "FHEt") <- "ct"
  attr(res, "FHEs") <- "FandV"
  res
}
# Again, see above for why this is here
`%*%.Rcpp_FandV_ct_mat` <- function(x, y) {
  if(class(y) == "Rcpp_FandV_ct_vec")
    y <- cbind(y)
  if(class(y) != "Rcpp_FandV_ct_mat")
    stop("requires cipher text matrix/vector arguments")
  if(x$ncol!=y$nrow) {
    stop("non-conformable arguments")
  }
  res <- x$matmulParallel(y)

  attr(res, "FHEt") <- "ctmat"
  attr(res, "FHEs") <- "FandV"
  res
}

loadFHE.Rcpp_FandV_ct <- function(file) {
  res <- load_FandV_ct(file, rlkLocker)
  attr(res, "FHEt") <- "ct"
  attr(res, "FHEs") <- "FandV"
  res
}

loadFHE.Rcpp_FandV_ct_vec <- function(file) {
  res <- load_FandV_ct_vec(file, rlkLocker)
  attr(res, "FHEt") <- "ctvec"
  attr(res, "FHEs") <- "FandV"
  res
}

loadFHE.Rcpp_FandV_ct_mat <- function(file) {
  res <- load_FandV_ct_mat(file, rlkLocker)
  attr(res, "FHEt") <- "ctmat"
  attr(res, "FHEs") <- "FandV"
  res
}

loadFHE.FandV_keys <- function(file) {
  res <- load_FandV_keys(file, rlkLocker)
  attr(res$pk, "FHEt") <- "pk"
  attr(res$pk, "FHEs") <- "FandV"
  attr(res$sk, "FHEt") <- "sk"
  attr(res$sk, "FHEs") <- "FandV"
  attr(res$rlk, "FHEt") <- "rlk"
  attr(res$rlk, "FHEs") <- "FandV"
  class(res) <- "FandV_keys"
  attr(res, "FHEt") <- "keys"
  attr(res, "FHEs") <- "FandV"
  res
}

# Dummies because:
#  1. Rcpp modules don't seem to work with S3method() in NAMESPACE
#  2. path.expand() in the saveFHE method doesn't overwrite the file argument
saveFHE.FandV_keys <- function(object, file) {
  saveFHE.FandV_keys2(object, path.expand(file))
}
saveFHE.Rcpp_FandV_ct <- function(object, file) {
  saveFHE.Rcpp_FandV_ct2(object, path.expand(file))
}
saveFHE.Rcpp_FandV_ct_vec <- function(object, file) {
  saveFHE.Rcpp_FandV_ct_vec2(object, path.expand(file))
}
saveFHE.Rcpp_FandV_ct_mat <- function(object, file) {
  saveFHE.Rcpp_FandV_ct_mat2(object, path.expand(file))
}
