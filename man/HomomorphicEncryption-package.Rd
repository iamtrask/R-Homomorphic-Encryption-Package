\name{HomomorphicEncryption-package}
\alias{HomomorphicEncryption-package}
\alias{HomomorphicEncryption}
\docType{package}
\title{
Homomorphic Encryption
}
\description{
This package provides easy to use implementations of some homomorphic encryption
schemes.  An encryption scheme is said to be homomorphic when certain functions
can be applied directly to the cipher text in such a way that decrypting the
result renders the same answer as if the function had been applied to the
unencrypted data.
}
\details{
\tabular{ll}{
Package: \tab HomomorphicEncryption\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2015-08-19\cr
License: \tab GPL-2\cr
}
}
\author{
Louis Aslett

Maintainer: Louis Aslett <aslett@stats.ox.ac.uk>
}
\references{
Aslett, L. J. M., Esperança, P. M. and Holmes, C. C. (2015), A review of homomorphic encryption and software tools for encrypted statistical machine learning. Technical report, University of Oxford.
}
\keyword{ package }
\examples{
# Generate cryptographic parameters
p <- pars("FandV")

# Create public/private keypair
keys <- keygen(p)

# Encrypt the values 2 and 3
ct1 <- enc(keys$pk, 2)
ct2 <- enc(keys$pk, 3)

# Homomorphically add the cipertexts together
ct <- ct1 + ct2

# Decrypt to 5, the result of applying + to plain messages
dec(keys$sk, ct)
}
