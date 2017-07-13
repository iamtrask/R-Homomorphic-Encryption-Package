# R-Homomorphic-Encryption-Package
This is a slightly modified fork of the R HomomorphicEncryption package originally from http://www.louisaslett.com/HomomorphicEncryption/

To use, clone the repository to your local disk. Open your "r" command line and run:

```
install.packages("/path/to/R-Homomorphic-Encryption-Package", repos = NULL, type="source")
library("HomomorphicEncryption")
p <- pars("FandV")
k <- keygen(p)
c1 <- enc(k$pk, c(42,34))
c2 <- enc(k$pk, c(7,5))
cres1 <- c1 + c2
cres2 <- c1 * c2
cres3 <- c1 %*% c2
dec(k$sk, cres1)
dec(k$sk, cres2)
dec(k$sk, cres3)
```
