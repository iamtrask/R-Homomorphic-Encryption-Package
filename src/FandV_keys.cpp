/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppParallel.h>
using namespace RcppParallel;

#include <limits.h>
#include <string>

#include "FandV_keys.h"
#include "FandV_ct.h"
#include "FandV_ct_vec.h"
#include "FandV_ct_mat.h"
#include "FandV.h"

#include <fmpz_polyxx.h>
using namespace flint;

//// Public keys ////
FandV_pk::FandV_pk(FandV_rlk_locker* rlkl, size_t rlki) : p(0, 0.0, 0, 1), rlkl(rlkl), rlki(rlki) { }

FandV_pk::FandV_pk(const FandV_pk& pk) : p(pk.p), rlkl(pk.rlkl), rlki(pk.rlki), p0(pk.p0), p1(pk.p1) { }

// Encrypt
void FandV_pk::enc(int m, FandV_ct& ct) const {
  ct.c0.realloc(p.Phi.length());
  ct.c1.realloc(p.Phi.length());
  
  fmpz_polyxx u, mP;
  u.realloc(p.Phi.length());
  mP.realloc(31);
  
  // Random numbers
  for(int i=0; i<p.Phi.length(); i++) {
    u.set_coeff(i, lround(R::rnorm(0.0,p.sigma))); // u
    ct.c0.set_coeff(i, lround(R::rnorm(0.0,p.sigma))); // e1
    ct.c1.set_coeff(i, lround(R::rnorm(0.0,p.sigma))); // e2
  }
  
  // Binary conversion of message
  int sign = 1;
  sign = copysign(sign, m);
  m = abs(m);
  for(int i=0; i<31; i++) {
    mP.set_coeff(i, (m&1)*sign);
    m >>= 1;
  }
  
  ct.c0 = ((p0*u)%p.Phi) + ct.c0 + p.Delta*mP;
  fmpz_polyxx_q(ct.c0, p.q);
  
  ct.c1 = ((p1*u)%p.Phi);
  fmpz_polyxx_q(ct.c1, p.q);
}

struct FandV_EncVec : public Worker {
  // Input values to encrypt & key
  const IntegerVector* input;
  const FandV_pk* pk;
  
  // Output vector of cipher texts
  std::vector<FandV_ct>* output;
  
  // Constructor
  FandV_EncVec(const FandV_pk* pk_, const IntegerVector* input_, std::vector<FandV_ct>* output_) { pk=pk_; input=input_; output=output_; }
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) {
      pk->enc((*input)[i], output->at(i));
    }
  }
};
void FandV_pk::encvec(IntegerVector m, FandV_ct_vec& ctvec) {
  FandV_ct ct(p, rlkl, rlki);
  ctvec.vec.resize(m.size(), ct);
  FandV_EncVec encEngine(this, &m, &(ctvec.vec));
  parallelFor(0, m.size(), encEngine);
}
void FandV_pk::encmat(IntegerVector m, int nrow, int ncol, FandV_ct_mat& ctmat) {
  FandV_ct ct(p, rlkl, rlki);
  ctmat.mat.resize(m.size(), ct);
  ctmat.nrow = nrow;
  ctmat.ncol = ncol;
  FandV_EncVec encEngine(this, &m, &(ctmat.mat));
  parallelFor(0, m.size(), encEngine);
}
  
void FandV_pk::show() {
  Rcout << "Fan and Vercauteren public key\n";
  Rcout << "( p\u2080 = ";
  printPoly(p0);
  Rcout << ",\np\u2081 = ";
  printPoly(p1);
  Rcout << " )\n";
}

// Save/load
void FandV_pk::save(FILE* fp) const {
  fprintf(fp, "=> FHE package object <=\nRcpp_FandV_pk\n");
  print(fp, p0);
  fprintf(fp, "\n");
  print(fp, p1);
  fprintf(fp, "\n");
}
FandV_pk::FandV_pk(FILE* fp, const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_) : p(p_), rlkl(rlkl_), rlki(rlki_) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE package object <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (PK)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_pk\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a public key\n";
    free(buf);
    return;
  }
  
  read(fp, p0);
  read(fp, p1);
  
  free(buf);
}


//// Private keys ////
FandV_sk::FandV_sk() { }

FandV_sk::FandV_sk(const FandV_sk& sk) : s(sk.s) { }

// Decrypt
std::vector<int> FandV_sk::decraw(const FandV_ct& ct) const {
  std::vector<int> result;
  fmpz_polyxx res, res2;
  fmpzxx tmp(1);
  
  res = ct.c0+((ct.c1*s)%ct.p.Phi);
  fmpz_polyxx_q(res, ct.p.q);
  res2 = (ct.p.t*res)%ct.p.q;
  res = (ct.p.t*res)/ct.p.q;
  for(int i=0; i<ct.p.Phi.length(); i++) {
    if(res2.get_coeff(i) > ct.p.q/2)
      res.set_coeff(i, res.get_coeff(i)+tmp);
  }
  fmpz_polyxx_q(res, ct.p.t);
  
  result.resize(res.length(), 0);
  for(int i=0; i<res.length(); i++) {
    result[i] = res.get_coeff(i).to<slong>();
  }
  
  return(result);
}
std::string FandV_sk::dec(const FandV_ct& ct) const {
  std::vector<int> res;
  fmpzxx tmp(1), m(0);
  
  res = decraw(ct);
  
  for(unsigned int i=0; i<res.size(); i++) {
    m += res[i]*tmp;
    tmp *= 2;
  }
  
  return(m.to_string());
}

void FandV_sk::show() {
  Rcout << "Fan and Vercauteren private key\n";
  Rcout << "s = ";
  printPoly(s);
  Rcout << "\n";
}

// Save/load
void FandV_sk::save(FILE* fp) const {
  fprintf(fp, "=> FHE package object <=\nRcpp_FandV_sk\n");
  print(fp, s);
  fprintf(fp, "\n");
}
FandV_sk::FandV_sk(FILE* fp) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE package object <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (SK)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_sk\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a secret key\n";
    free(buf);
    return;
  }
  
  read(fp, s);
  free(buf);
}


//// Relinearisation keys ////
FandV_rlk::FandV_rlk() { }

FandV_rlk::FandV_rlk(const FandV_rlk& rlk) : rlk00(rlk.rlk00), rlk01(rlk.rlk01), rlk10(rlk.rlk10), rlk11(rlk.rlk11) { }

void FandV_rlk::show() {
  Rcout << "Fan and Vercauteren relinearisation key\n";
  Rcout << "( rlk\u2080\u2080 = ";
  printPoly(rlk00);
  Rcout << ",\nrlk\u2080\u2081 = ";
  printPoly(rlk01);
  Rcout << ",\nrlk\u2081\u2080 = ";
  printPoly(rlk10);
  Rcout << ",\nrlk\u2081\u2081 = ";
  printPoly(rlk11);
  Rcout << " )\n";
}

// Save/load
void FandV_rlk::save(FILE* fp) const {
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_rlk\n");
  print(fp, rlk00);
  fprintf(fp, "\n");
  print(fp, rlk01);
  fprintf(fp, "\n");
  print(fp, rlk10);
  fprintf(fp, "\n");
  print(fp, rlk11);
  fprintf(fp, "\n");
}
FandV_rlk::FandV_rlk(FILE* fp) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (RLK)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_rlk\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a relinearisation key\n";
    free(buf);
    return;
  }
  
  read(fp, rlk00);
  read(fp, rlk01);
  read(fp, rlk10);
  read(fp, rlk11);
  
  free(buf);
}

FandV_rlk_locker::FandV_rlk_locker() { }

int FandV_rlk_locker::add(const FandV_rlk &rlk) {
  for(unsigned int i=0; i<x.size(); i++) {
    if(x[i].rlk00 == rlk.rlk00 && x[i].rlk01 == rlk.rlk01 && x[i].rlk10 == rlk.rlk10 && x[i].rlk11 == rlk.rlk11) {
      return(i);
    }
  }
  x.push_back(rlk);
  return(x.size()-1);
}

void FandV_rlk_locker::show() const {
  Rcout << "Locker contains " << x.size() << " Fan and Vercauteren relinearisation keys\n";
}
