/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/
//#define RCPP_DEBUG_LEVEL 1
#include <Rcpp.h>
using namespace Rcpp;

#include <fmpzxx.h>
#include <fmpz_polyxx.h>
using namespace flint;

#include <vector>
#include <stdio.h>

#include "FandV_par.h"
#include "FandV_keys.h"
#include "FandV_ct.h"
#include "FandV_ct_vec.h"
#include "FandV_ct_mat.h"

// More detailed info on memory usage.  Rcpp modules exist outside R's direct
// control, so gc() useless for finding out memory usage.
#include "../config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
void HEmem() {
  #ifdef HAVE_MALLOC_STATS
    malloc_stats();
  #else
    Rcout << "malloc_stats() not available on this machine\n";
  #endif
}

// Do centred modulo q reduction of all coefficients of polynomial p ... [p]_q
void fmpz_polyxx_q(fmpz_polyxx& p, fmpzxx q) {
  fmpzxx tmp, qo2(q/2);
  for(int i=0; i<p.length(); i++) {
    tmp = p.get_coeff(i)%q;
    if(tmp > qo2)
      tmp -= q;
    p.set_coeff(i, tmp);
  }
}

void fmpz_rand(fmpzxx &p, unsigned int bits) { // Random number from 0 to 2^bits-1
  RNGScope scope;
  p = 0;
  for(unsigned int i=0; i<bits/32; i++) {
    p = (p << 32) + ((unsigned int) R::runif(0.0, 4294967295.0));
  }
  if(bits%32 > 0) {
    p = (p << (bits%32)) + ((unsigned int) R::runif(0.0, pow(2.0, (double) (bits%32))-1.0));
  }
}

void printPoly(const fmpz_polyxx& p) {
  static const char * const super[] = {"\xe2\x81\xb0", "\xc2\xb9", "\xc2\xb2",
    "\xc2\xb3", "\xe2\x81\xb4", "\xe2\x81\xb5", "\xe2\x81\xb6",
    "\xe2\x81\xb7", "\xe2\x81\xb8", "\xe2\x81\xb9"};
  
  bool firstDone=false;
  for(int i=p.length(); --i>=0; ) {
    if(p.get_coeff(i) == fmpzxx(0)) continue;
    if(p.get_coeff(i)>0 && firstDone && i<p.length()) Rcout << "+";
    if(p.get_coeff(i) == fmpzxx(-1) && i>0)
      Rcout << "-";
    else if(p.get_coeff(i) != fmpzxx(1) || i==0)
      Rcout << p.get_coeff(i);
    firstDone=true;
    if(i > 0)
      Rcout << "x";
    if(i > 1) {
      int j, max;
      j = i;
      max = 1;
      while(j > 9) {
        max *= 10;
        j /= 10;
      }
      j = i;
      for(int k=max; k>0; k/=10)
        Rcout << super[(j/k)%10];
    }
  }
}

void save_FandV_ct(const FandV_ct& ct, const std::string& file) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "w");
  if(fp == NULL) {
    perror("Error");
  }
  
  // header
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_ct\n");
  // pars
  ct.p.save(fp);
  // rlk
  (ct.rlkl->x[ct.rlki]).save(fp);
  // ct content
  ct.save(fp);
  
  fclose(fp);
}
FandV_ct load_FandV_ct(const std::string& file, FandV_rlk_locker* rlkl) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "r");
  if(fp == NULL) {
    perror("Error");
  }

  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (CT (a))\n";
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_ct\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a single ciphertext object\n";
  }

  // pars
  FandV_par p(fp);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // rlk
  FandV_rlk rlk(fp);
  int rlki = rlkl->add(rlk);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // ct content
  FandV_ct ct(fp, p, rlkl, rlki);
  
  fclose(fp);
  free(buf);
  
  return(ct);
}

void save_FandV_ct_vec(const FandV_ct_vec& ct_vec, const std::string& file) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "w");
  if(fp == NULL) {
    perror("Error");
  }
  
  // header
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_ct_vec\n");
  // pars
  ct_vec.vec[1].p.save(fp);
  // rlk
  (ct_vec.vec[1].rlkl->x[ct_vec.vec[1].rlki]).save(fp);
  // ct content
  ct_vec.save(fp);
  
  fclose(fp);
}
FandV_ct_vec load_FandV_ct_vec(const std::string& file, FandV_rlk_locker* rlkl) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "r");
  if(fp == NULL) {
    perror("Error");
  }

  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (VEC (a))\n";
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_ct_vec\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a ciphertext vector object\n";
  }

  // pars
  FandV_par p(fp);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // rlk
  FandV_rlk rlk(fp);
  int rlki = rlkl->add(rlk);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // ct_vec content
  FandV_ct_vec ct_vec(fp, p, rlkl, rlki);
  
  fclose(fp);
  free(buf);
  
  return(ct_vec);
}

void save_FandV_ct_mat(const FandV_ct_mat& ct_mat, const std::string& file) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "w");
  if(fp == NULL) {
    perror("Error");
  }
  
  // header
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_ct_mat\n");
  // pars
  ct_mat.mat[1].p.save(fp);
  // rlk
  (ct_mat.mat[1].rlkl->x[ct_mat.mat[1].rlki]).save(fp);
  // ct content
  ct_mat.save(fp);
  
  fclose(fp);
}
FandV_ct_mat load_FandV_ct_mat(const std::string& file, FandV_rlk_locker* rlkl) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "r");
  if(fp == NULL) {
    perror("Error");
  }

  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (MAT (a))\n";
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_ct_mat\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a ciphertext matrix object\n";
  }

  // pars
  FandV_par p(fp);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // rlk
  FandV_rlk rlk(fp);
  int rlki = rlkl->add(rlk);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // ct_vec content
  FandV_ct_mat ct_mat(fp, p, rlkl, rlki);
  
  fclose(fp);
  free(buf);
  
  return(ct_mat);
}

void save_FandV_keys(const List& keys, const std::string& file) {
  const char *file_c = file.c_str();
  
  FILE *fp = fopen(file_c, "w");
  if(fp == NULL) {
    perror("Error");
  }
  
  fprintf(fp, "=> FHE pkg obj <=\nFandV_keys\n");
  
  FandV_rlk rlk = keys["rlk"];
  FandV_pk pk = keys["pk"];
  FandV_sk sk = keys["sk"];
  
  // pars + rlk
  pk.p.save(fp);
  rlk.save(fp);
  // pk
  pk.save(fp);
  // sk
  sk.save(fp);
  
  fclose(fp);
}
List load_FandV_keys(const std::string& file, FandV_rlk_locker* rlkl) {
  const char *file_c = file.c_str();
  List keys;
  
  FILE *fp = fopen(file_c, "r");
  if(fp == NULL) {
    perror("Error");
  }
  
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (KEYS)\n";
    free(buf);
    return(keys);
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("FandV_keys\n", buf, len) != 0) {
    Rcout << "Error: file does not contain key objects\n";
    free(buf);
    return(keys);
  }
  
  // pars
  FandV_par p(fp);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // rlk
  FandV_rlk rlk(fp);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  int rlki = rlkl->add(rlk);
  // pk
  FandV_pk pk(fp, p, rlkl, rlki);
  len = getline(&buf, &bufn, fp); // Advance past the new line
  // sk
  FandV_sk sk(fp);
  
  keys["sk"] = sk;
  keys["pk"] = pk;
  keys["rlk"] = rlk;
  
  fclose(fp);
  
  free(buf);
  return(keys);
}


RCPP_EXPOSED_CLASS(FandV_par)
RCPP_EXPOSED_CLASS(FandV_pk)
RCPP_EXPOSED_CLASS(FandV_sk)
RCPP_EXPOSED_CLASS(FandV_rlk)
RCPP_EXPOSED_CLASS(FandV_rlk_locker)
RCPP_EXPOSED_CLASS(FandV_ct)
RCPP_EXPOSED_CLASS(FandV_ct_vec)
RCPP_EXPOSED_CLASS(FandV_ct_mat)

RCPP_MODULE(FandV) {
  class_<FandV_par>("FandV_par")
    .constructor<int, double, int, int, int, int>()
    .method("keygen", &FandV_par::keygen)
    .method("show", &FandV_par::show)
    .method("show_no_t", &FandV_par::show_no_t)
    .method("show_t", &FandV_par::show_t)
    .method("get_t", &FandV_par::get_t)
  ;
  
  class_<FandV_pk>("FandV_pk")
    .constructor<FandV_rlk_locker*, size_t>()
    .field("p", &FandV_pk::p)
    .field("rlki", &FandV_pk::rlki)
    .method("enc", &FandV_pk::enc)
    .method("encvec", &FandV_pk::encvec)
    .method("encmat", &FandV_pk::encmat)
    .method("show", &FandV_pk::show)
  ;

  class_<FandV_sk>("FandV_sk")
    .constructor()
    .method("decraw", &FandV_sk::decraw)
    .method("dec", &FandV_sk::dec)
    .method("show", &FandV_sk::show)
  ;
  
  class_<FandV_rlk>("FandV_rlk")
    .constructor()
    .method("show", &FandV_rlk::show)
  ;
  
  class_<FandV_rlk_locker>("FandV_rlk_locker")
    .constructor()
    .method("add", &FandV_rlk_locker::add)
    .method("show", &FandV_rlk_locker::show)
  ;
  
  class_<FandV_ct>("FandV_ct")
    .constructor<FandV_par,FandV_rlk_locker*,size_t>()
    .field("p", &FandV_ct::p)
    .field("rlki", &FandV_ct::rlki)
    .field("depth", &FandV_ct::depth)
    .method("add", &FandV_ct::add)
    .method("sub", &FandV_ct::sub)
    .method("mul", &FandV_ct::mul)
    .method("show", &FandV_ct::show)
  ;
  
  class_<FandV_ct_vec>("FandV_ct_vec")
    .constructor()
    .method("add", &FandV_ct_vec::add)
    .method("addct", &FandV_ct_vec::addct)
    .method("sub", &FandV_ct_vec::sub)
    .method("subct", &FandV_ct_vec::subct)
    .method("get", &FandV_ct_vec::get)
    .method("mulParallel", &FandV_ct_vec::mulParallel)
    .method("mulSerial", &FandV_ct_vec::mulSerial)
    .method("mulctParallel", &FandV_ct_vec::mulctParallel)
    .method("mulctSerial", &FandV_ct_vec::mulctSerial)
    .method("sumParallel", &FandV_ct_vec::sumParallel)
    .method("sumSerial", &FandV_ct_vec::sumSerial)
    .method("prodParallel", &FandV_ct_vec::prodParallel)
    .method("prodSerial", &FandV_ct_vec::prodSerial)
    .method("innerprod", &FandV_ct_vec::innerprod)
    .method("push", &FandV_ct_vec::push)
    .method("pushvec", &FandV_ct_vec::pushvec)
    .method("set", &FandV_ct_vec::set)
    .method("show", &FandV_ct_vec::show)
    .method("size", &FandV_ct_vec::size)
    .method("subset", &FandV_ct_vec::subset)
    .method("without", &FandV_ct_vec::without)
  ;
  
  class_<FandV_ct_mat>("FandV_ct_mat")
    .constructor()
    .field("nrow", &FandV_ct_mat::nrow)
    .field("ncol", &FandV_ct_mat::ncol)
    .method("size", &FandV_ct_mat::size)
    .method("get", &FandV_ct_mat::get)
    .method("subset", &FandV_ct_mat::subset)
    .method("subsetV", &FandV_ct_mat::subsetV)
    .method("t", &FandV_ct_mat::t)
    .method("set", &FandV_ct_mat::set)
    .method("setelt", &FandV_ct_mat::setelt)
    .method("setmatrix", &FandV_ct_mat::setmatrix)
    .method("reset", &FandV_ct_mat::reset)
    .method("show", &FandV_ct_mat::show)
    .method("add", &FandV_ct_mat::add)
    .method("mul", &FandV_ct_mat::mul)
    .method("addct", &FandV_ct_mat::addct)
    .method("mulctParallel", &FandV_ct_mat::mulctParallel)
    .method("mulctSerial", &FandV_ct_mat::mulctSerial)
    .method("mulctvecParallel", &FandV_ct_mat::mulctvecParallel)
    .method("mulctvecSerial", &FandV_ct_mat::mulctvecSerial)
    .method("matmulParallel", &FandV_ct_mat::matmulParallel)
    .method("matmulSerial", &FandV_ct_mat::matmulSerial)
    .method("TmatmulParallel", &FandV_ct_mat::TmatmulParallel)
    .method("matmulTParallel", &FandV_ct_mat::matmulTParallel)
    .method("rowSumsParallel", &FandV_ct_mat::rowSumsParallel)
    .method("rowSumsSerial", &FandV_ct_mat::rowSumsSerial)
    .method("colSumsParallel", &FandV_ct_mat::colSumsParallel)
    .method("colSumsSerial", &FandV_ct_mat::colSumsSerial)
  ;
  
  function("saveFHE.FandV_keys2", &save_FandV_keys);
  function("load_FandV_keys", &load_FandV_keys);
  function("saveFHE.Rcpp_FandV_ct2", &save_FandV_ct);
  function("load_FandV_ct", &load_FandV_ct);
  function("saveFHE.Rcpp_FandV_ct_vec2", &save_FandV_ct_vec);
  function("load_FandV_ct_vec", &load_FandV_ct_vec);
  function("saveFHE.Rcpp_FandV_ct_mat2", &save_FandV_ct_mat);
  function("load_FandV_ct_mat", &load_FandV_ct_mat);
  function("HEmem", &HEmem);
}
