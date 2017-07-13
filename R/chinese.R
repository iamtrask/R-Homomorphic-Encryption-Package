# Supporting functions for private namespace in the package.  We need bigint
# versions of these, which are adapted from the numbers package
chinese <- function(a, m) {
  stopifnot(is.bigz(a), is.bigz(m))
  n <- length(m)
  if(length(a) != n) 
    stop("arguments 'a' an 'm' must be vectors of equal length.")
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      if(gcd(m[i], m[j]) > 1) {
        stop("elements in argument 'm' must be pairwise coprime.")
      }
    }
  }
  M <- prod(m)
  x <- 0
  for(i in 1:n) {
    Mmi <- prod(m[-i])
    mmi <- modinv(Mmi, m[i])
    x <- x + a[i] * Mmi * mmi
  }
  return(x%%M)
}

modinv <- function(n, m) {
  stopifnot(is.bigz(n), is.bigz(m))
  v <- extGCD(n, m)
  if(v[1] == as.bigz(0) || v[1] > as.bigz(1))
    return(NA)
  if(v[2] >= as.bigz(0))
    v[2]
  else v[2] + m
}

extGCD <- function(a, b) {
  stopifnot(is.bigz(a), length(a) == 1,
            is.bigz(b), length(b) == 1)
  sign_ab <- sign(c(a, b))
  A <- matrix(c(abs(c(a, b)), 1, 0, 0, 1), nrow = 2, ncol = 3)
  while(A[1, 1] * A[2, 1] != as.bigz(0)) {
    if(A[1, 1] > A[2, 1]) {
      m <- A[1, 1]%/%A[2, 1]
      A[1, ] <- A[1, ] - m * A[2, ]
    } else {
      m <- A[2, 1]%/%A[1, 1]
      A[2, ] <- A[2, ] - m * A[1, ]
    }
  }
  if(A[1, 1] == as.bigz(0)) 
    g <- c(A[2, ]) # bigz ignores drop=TRUE
  else
    g <- c(A[1, ]) # bigz ignores drop=TRUE
  g[2:3] <- sign_ab * g[2:3]
  return(g)
}
