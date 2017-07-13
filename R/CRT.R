# schemeName <- function(name) {
#   switch(name, "FandV_CRT"="Fan and Vercauteren")
# }
# 
# evalqOnLoad({
#   setClass("CRT", slots=c(ct="list", ctvalid="logical"))
#   
#   setMethod("show", signature(object="CRT"), function(object) {
#     if(attr(object, "FHEt")=="ct")
#       cat(schemeName(attr(object, "FHEs")), " cipher text with Chinese Remainder Theorem message space modulus extension (", sum(object@ctvalid)-1, " comparisons remaining)\n", sep="")
#     else if(attr(object, "FHEt") == "ctvec")
#       cat("Vector of", length(object@ct[[1]]), schemeName(attr(object, "FHEs")), "cipher texts with Chinese Remainder Theorem message space modulus extension\n")
#     else if(attr(object, "FHEt") == "ctmat")
#       cat("Matrix of", dim(object@ct[[1]])[1],"x", dim(object@ct[[1]])[2], schemeName(attr(object, "FHEs")), "cipher texts with Chinese Remainder Theorem message space modulus extension\n")
#   })
#   
#   setMethod("Ops", signature(e1="CRT", e2="CRT"), function(e1, e2) {
#     crt <- new("CRT", ct=list())
#     
#     for(i in 1:length(e1@ct)) {
#       crt@ct[[i]] <- callGeneric(e1@ct[[i]], e2@ct[[i]])
#     }
#     
#     attr(crt, "FHEt") <- attr(crt@ct[[1]], "FHEt")
#     attr(crt, "FHEs") <- "FandV_CRT"
#     return(crt)
#   })
# })
