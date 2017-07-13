#' Encryption parameter selection
#' 
#' Use this function to help choose encryption parameters for your scheme.
#' 
#' This function provides assistance with selecting parameter values for the
#' encryption scheme you wish to use.
#' 
#' Currently only the scheme of Fan and Vercauteren (\code{"FandV"}) is implemented.
#' 
#' For \code{"FandV"} you may specify:
#' \describe{
#'   \item{\code{lambda}}{the security level required (bits), default is 80;}
#'   \item{\code{max}}{the largest absolute value you will need to store encrypted,
#'    default is 1000;}
#'   \item{\code{L}}{the deepest multiplication depth you need to be able to 
#'   evaluate encrypted, default is 4.}
#' }
#' 
#' The security level in bits relates to the approximate number of operations
#' which would be required to break the cipher text, i.e. \eqn{O(2^\lambda)}
#' operations.
#' 
#' Ensuring the depth of multiplication operations you require is based on an
#' overwhelming probability of being able to correctly decrypt after so many
#' operations by roughly bounding the possible noise growth.  In most situations
#' it is very conservative and additional multiplications are possible, so this
#' should be considered only a lower bound on the number of multiplies required.  If 
#' computational performance is a concern this should be treated only as a good starting 
#' point from which the relevant parameters can be eased off, testing more thoroughly
#' by simulation.
#' 
#' For the scheme of Fan and Vercauteren (\code{"FandV"}): the security level is
#' based on the analysis of Lindner and Peikert (2011) which is known to be quite
#' conservative; and the multiplicative depth bound is based on the analysis of
#' Lepoint and Naehrig (2014).  If the multiplicative depth lower bound is too
#' conservative, then for computational performance reduce the qpow parameter from
#' this recommended starting point.
#' 
#' @param scheme the scheme for which to get parameter help.  Currently
#' only Fan and Vercauteren's scheme is supported by specifying \code{"FandV"}.
#' 
#' @param ... pass the specific options for the chosen scheme as named arguments.
#' See the details section for options for encryption schemes currently implemented.
#' 
#' @references
#' Fan, J., & Vercauteren, F. (2012). Somewhat Practical Fully Homomorphic
#' Encryption. IACR ePrint. Retrieved from \url{https://eprint.iacr.org/2012/144}
#' 
#' Lepoint, T., & Naehrig, M. (2014). A comparison of the homomorphic encryption 
#' schemes FV and YASHE. In \emph{Progress in Cryptology–AFRICACRYPT 2014} 
#' (pp. 318-335).
#' 
#' Lindner, R., & Peikert, C. (2011). Better key sizes (and attacks) for LWE-based
#' encryption. In \emph{Topics in Cryptology–CT-RSA 2011} (pp. 319-339).
#' 
#' @seealso
#' \code{\link{pars}} to manually choose parameters.
#' 
#' \code{\link{keygen}} to generate encryption keys using these parameters.
#' 
#' @examples
#' # Want 128-bit security with 6 multiplies deep
#' # to evaluate 7 factorial
#' p <- parsHelp("FandV", lambda=128, L=6)
#' keys <- keygen(p)
#' ct <- enc(keys$pk, 1:7)
#' dec(keys$sk, prod(ct))
#' factorial(7)
#' 
#' # But notice this is quite conservative
#' # This will usually give the right answer too
#' ct <- enc(keys$pk, 1:8)
#' dec(keys$sk, prod(ct))
#' factorial(8)
#' 
#' @author Louis Aslett
parsHelp <- function(scheme, ...) {
  args <- list(...)
  if(scheme=="FandV") {
    lambda <- 80
    max <- 1000
    L <- 4
    if("lambda" %in% names(args)) {
      if(lambda<32) stop("Less than 32 bits of security is very weak.")
      lambda <- args[["lambda"]]
    }
    if("max" %in% names(args)) {
      if(max<1) stop("Absolute value of maximum value to store cannot be less than 1.")
      max <- args[["max"]]
    }
    if("L" %in% names(args)) {
      if(L<0) stop("Cannot specify negative multiplicative depth.")
      L <- args[["L"]]
    }
    
    # Sensible defaults for the unspecified parameters
    sigma <- 16
    B_err <- as.bigz(3*sigma)
    
    # Start parameter hunt
    dpow <- 8
    
    maxqpow <- FandV_maxqLP11(2^dpow, lambda, sigma)
    while(!FandV_testDepth(L, d=as.bigz(2)^dpow, q=as.bigz(2)^maxqpow, t=max*2, B_err=B_err)) {
      dpow <- dpow+1
      maxqpow <- floor(FandV_maxqLP11(2^dpow, lambda, sigma))
      if(maxqpow%%2!=0) { # qpow must be even
        maxqpow <- maxqpow-1
      }
    }
    # We now have the value of dpow, but maxqpow might be bigger than necessary
    # so we start to pull it down to see where the tipping point is
    while(FandV_testDepth(L, d=as.bigz(2)^dpow, q=as.bigz(2)^maxqpow, t=max*2, B_err=B_err)) {
      maxqpow <- maxqpow - 2 # qpow must be even
    }
    qpow <- maxqpow + 2
    # We now have parameters!
    return(pars("FandV", d=2^dpow, sigma=sigma, qpow=qpow, t=max*2))
  }
  stop("The scheme ", scheme, " is not recognised.  Currently only 'FandV' is implemented.")
}

# Returns true if parameter set can support multiplicative depth L with
# overwhelming probability
FandV_testDepth <- function(L, d=as.bigz(2)^11, q=as.bigz(2)^80, t=as.bigz(2)^10, B_key=as.bigz(1), B_err=as.bigz(8), w=as.bigz(2^32)) {
  C1 <- function(delta, t, B_key) {
    #(1+4*(delta*B_key)^{-1})*delta*t*B_key
    delta^2*t*B_key+4*delta*t
  }
  
  C2 <- function(delta, t, B_key, B_err, lwq, w) {
    delta^2*B_key*(B_key+t^2) + delta*lwq*w*B_err
  }
  
  V <- function(delta, B_key, B_err) {
    B_err*(1+2*delta*B_key)
  }
  
  rhs <- function(q, t) {
    2^{log2(q)-log2(t)-1}
  }
  
  lwq <- floor(log(q, w)) + 1
  delta <- as.bigz(sqrt(as.integer(d))) # possibly quite weak upper bound based on Cauchy Schwartz
  LHS <- C1(delta, t, B_key)^L * V(delta, B_key, B_err) + L * C1(delta, t, B_key)^{L-1} * C2(delta, t, B_key, B_err, lwq, w)
  I(rhs(q, t)-LHS>0)
}

# Find maximum bit length of q to achieve required security using Lepoint and
# Naehrig 2014
# Not used at present
#   epsilon = adversary's advantage
FandV_maxqLN14 <- function(d, lambda, sigma, epsilon=2^{-64}) {
  alpha <- function(epsilon) {
    sqrt(-log10(epsilon)/pi)
  }
  # NB This is not the minimal root Hermite factor from their paper -- haven't
  #    worked out how to compute their table values for all possible lambda and m
  #    yet.  This is the value from F&V's paper
  gamma <- function(lambda) {
    2^(1.8/(lambda+110))
    #    1.00799
  }
  optimise(function(m) { (m^2 * log2(gamma(lambda)) + m * log2(sigma/alpha(epsilon)))/(m-d) }, interval=c(d+1, d^2))
}

# Find maximum bit length of q to achieve required security using Lindner & Peikert
#   epsilon = adversary's advantage
FandV_maxqLP11 <- function(d, lambda, sigma, epsilon=2^{-64}) {
  l2delta <- 1.8/(lambda+110)
  l2alpha <- log2(sqrt(log(1/epsilon)/pi))
  uniroot(function(x) { 2*sqrt(d*x*l2delta)-l2alpha-x+log2(sigma) }, interval=c(1,10000000))$root
}

