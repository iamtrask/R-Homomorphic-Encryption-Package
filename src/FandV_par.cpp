/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#include <Rcpp.h>
using namespace Rcpp;

#include <arith.h>

#include "FandV_par.h"
#include "FandV.h"
#include "FandV_keys.h"

// Construct from parameters
FandV_par::FandV_par(int d_, double sigma_, int qpow_, int t_, int lambda_, int L_) : sigma(sigma_), qpow(qpow_), q(1), t(t_), T(1), lambda(lambda_), L(L_) {
  arith_cyclotomic_polynomial(Phi._data().inner, 2*d_); // Phi is 2d-th Cyclotomic polynomial
  
  q = q << qpow; // q=2^qpow
  Delta = q/t;
  T = T << (qpow/2);
}

// Copy constructor
FandV_par::FandV_par(const FandV_par& par) : sigma(par.sigma), qpow(par.qpow), q(par.q), t(par.t), T(par.T), Delta(par.Delta), Phi(par.Phi), lambda(par.lambda), L(par.L) { }

// Swap function
void FandV_par::swap(FandV_par& a, FandV_par& b) {
  std::swap(a.sigma, b.sigma);
  std::swap(a.qpow, b.qpow);
  std::swap(a.q, b.q);
  std::swap(a.t, b.t);
  std::swap(a.T, b.T);
  std::swap(a.Delta, b.Delta);
  std::swap(a.Phi, b.Phi);
  std::swap(a.lambda, b.lambda);
  std::swap(a.L, b.L);
}

// Assignment (copy-and-swap idiom)
FandV_par& FandV_par::operator=(FandV_par par) {
  swap(*this, par);
  return(*this);
}

// Print
void FandV_par::show() {
  Rcout << "Fan and Vercauteren parameters\n";
  Rcout << "\u03d5 = ";
  printPoly(Phi);
  Rcout << "\nq = " << q << " (" << qpow << "-bit integer)\nt = " << t << "\n\u0394 = " << Delta << "\n\u03c3 = " << sigma << "\nSecurity level \u2248 " << lambda << "-bits\nSupports multiplicative depth of " << L << " with overwhelming probability (i.e. lower bound, likely more possible)\n";
}
void FandV_par::show_no_t() {
  Rcout << "\u03d5 = ";
  printPoly(Phi);
  Rcout << "\nq = " << q << " (" << qpow << "-bit integer)\n\u0394 = " << Delta << "\n\u03c3 = " << sigma << "\n";
}
void FandV_par::show_t() {
  Rcout << t;
}
std::string FandV_par::get_t() {
  return(t.to_string());
}

// Keygen
void FandV_par::keygen(FandV_pk& pk, FandV_sk& sk, FandV_rlk& rlk) {
  // WARNING: according to flint.h, flint_randinit() uses a fixed seed.
  //   https://github.com/wbhart/flint2/issues/93
  //frandxx fr;
  // WARNING: randtest() does strange things ... get inflated number of zeros etc
  //sk.sk = fmpz_polyxx::randtest_unsigned(fr, p.Phi().length(), (mp_bitcnt_t) 1);
  //pk.pk1 = fmpz_polyxx::randtest(fr, p.Phi().length(), (mp_bitcnt_t) (p.qpow()-1));
  // Use GMP and conversion instead
  // WARNING: fmpzxx and fmp_polyxx can't convert from GMP, so have to use C interface
  
  RNGScope scope;
  
  // Public/private keys
  pk.p = *this;
  
  fmpz_polyxx e;
  
  fmpzxx tmp, qo2p1(1);
  qo2p1 = (qo2p1 << (pk.p.qpow-1)) + fmpzxx(1);
  
  // Size up the polynomials
  sk.s.realloc(pk.p.Phi.length());
  pk.p0.realloc(pk.p.Phi.length());
  pk.p1.realloc(pk.p.Phi.length());
  e.realloc(pk.p.Phi.length());
  
  // Generate random parts
  for(unsigned int i=0; i<pk.p.Phi.length(); i++) {
    // s
    sk.s.set_coeff(i, lround(R::runif(0.0,1.0)));
    
    // a
    fmpz_rand(tmp, pk.p.qpow); // tmp \in (0, 2^q-1)
    tmp -= qo2p1; // tmp - 2^{q-1} + 1 \in (-2^{q-1}, 2^{q-1}]
    pk.p0.set_coeff(i, tmp);
    pk.p1.set_coeff(i, tmp);
    
    // e
    e.set_coeff(i, lround(R::rnorm(0.0,pk.p.sigma)));
  }
  // -(a.s+e) ...
  pk.p0 = -( ((pk.p0*sk.s)%pk.p.Phi) + e );
  // ... mod q
  fmpz_polyxx_q(pk.p0, pk.p.q);
  
  // Relin key
  for(unsigned int i=0; i<pk.p.Phi.length(); i++) {
    // a0
    fmpz_rand(tmp, pk.p.qpow); // tmp \in (0, 2^q-1)
    tmp -= qo2p1; // tmp - 2^{q-1} + 1 \in (-2^{q-1}, 2^{q-1}]
    rlk.rlk01.set_coeff(i, tmp);
    // a1
    fmpz_rand(tmp, pk.p.qpow); // tmp \in (0, 2^q-1)
    tmp -= qo2p1; // tmp - 2^{q-1} + 1 \in (-2^{q-1}, 2^{q-1}]
    rlk.rlk11.set_coeff(i, tmp);
    
    // e
    rlk.rlk00.set_coeff(i, lround(R::rnorm(0.0,pk.p.sigma)));
    rlk.rlk10.set_coeff(i, lround(R::rnorm(0.0,pk.p.sigma)));
  }
  // e var will now hold s^2
  e = ((sk.s*sk.s)%pk.p.Phi);
  rlk.rlk00 = -( ((rlk.rlk01*sk.s)%pk.p.Phi) + rlk.rlk00 ) + e;
  fmpz_polyxx_q(rlk.rlk00, pk.p.q);
  rlk.rlk10 = -( ((rlk.rlk11*sk.s)%pk.p.Phi) + rlk.rlk10 ) + T*e;
  fmpz_polyxx_q(rlk.rlk00, pk.p.q);
  
  // Make sure public key holds a copy of rlk so it can be passed onto ciphertexts
  pk.rlki = pk.rlkl->add(rlk);
}

// Save/load
void FandV_par::save(FILE* fp) const {
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_par\n");
  
  // sigma, qpow
  fprintf(fp, "%f:%d\n", sigma, qpow);
  
  // q
  print(fp, q);
  fprintf(fp, "\n");
  // t
  print(fp, t);
  fprintf(fp, "\n");
  // T
  print(fp, T);
  fprintf(fp, "\n");
  // Delta
  print(fp, Delta);
  fprintf(fp, "\n");
  // Phi
  print(fp, Phi);
  fprintf(fp, "\n");
}
FandV_par::FandV_par(FILE* fp) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (PAR)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_par\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a parameters object\n";
    free(buf);
    return;
  }
  
  // sigma, qpow
  len = fscanf(fp, "%lf:%d\n", &sigma, &qpow);
  // q
  read(fp, q);
  // t
  read(fp, t);
  // T
  read(fp, T);
  // Delta
  read(fp, Delta);
  // Phi
  read(fp, Phi);
  
  free(buf);
}
