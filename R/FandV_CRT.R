# #' Fan and Vercauteren encryption scheme with Chinese Remainder Theorem extension
# #' 
# #' The Fan and Vercauteren scheme is implemented in this package, together with a
# #' seamless implementation of modulus extension via the Chinese Remainder Theorem.
# #' 
# #' Description of the scheme.
# #' 
# #' @name FandV_CRT
# #' 
# #' @usage
# #' p <- pars("FandV_CRT")
# #' 
# #' @examples
# #' # Benchmark the performance of the scheme
# #' #library(microbenchmark)
# #' #p <- pars("FandV_CRT")
# #' #microbenchmark({ keys <- keygen(p) }, unit="ms")
# #' #microbenchmark({ ct1 <- enc(keys$pk, 2) }, unit="ms")
# #' #ct2 <- enc(keys$pk, 3)
# #' #microbenchmark({ ct1 + ct2 }, unit="ms")
# #' #microbenchmark({ ct1 * ct2 }, unit="ms")
# #' #microbenchmark({ dec(keys$sk, ct1) }, unit="ms")
# #' 
# NULL
# 
# # PRINTING COMMANDS
# # Parameters
# print.FandV_CRT <- function(x, ...) {
#   p <- x
#   cat("Fan and Vercauteren with Chinese Remainder Theorem message space modulus extension\n")
#   p[[1]]$show_no_t()
#   cat("with coprime moduli:\nt = {")
#   for(i in 1:(length(p)-1)) {
#     p[[i]]$show_t()
#     cat(", ")
#   }
#   p[[length(p)]]$show_t()
#   cat("}\n")
# }
# 
# # Keys
# print.FandV_CRT_keys <- function(x, ...) {
#   cat("Set of",length(x$pk),"keys for Fan and Vercauteren with Chinese Remainder Theorem message space modulus extension\n")
# }
# print.FandV_CRT_pk <- function(x, ...) {
#   cat("Set of",length(x),"public keys for Fan and Vercauteren with Chinese Remainder Theorem message space modulus extension\n")
# }
# print.FandV_CRT_sk <- function(x, ...) {
#   cat("Set of",length(x),"secret keys for Fan and Vercauteren with Chinese Remainder Theorem message space modulus extension\n")
# }
# print.FandV_CRT_rlk <- function(x, ...) {
#   cat("Set of",length(x),"relinearisation keys for Fan and Vercauteren with Chinese Remainder Theorem message space modulus extension\n")
# }
# 
# # Cipher texts
# 
