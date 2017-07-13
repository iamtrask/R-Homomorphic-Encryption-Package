/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#include <Rcpp.h>
using namespace Rcpp;

#include <RcppParallel.h>
using namespace RcppParallel;

#include <iostream>
#include <math.h>

#include "FandV_ct.h"
#include "FandV_ct_vec.h"
#include "FandV.h"

// Construct from parameters
FandV_ct_vec::FandV_ct_vec() { }
FandV_ct_vec::FandV_ct_vec(const std::vector<FandV_ct> v) : vec(v) { }

// Copy constructor
FandV_ct_vec::FandV_ct_vec(const FandV_ct_vec& ct_vec) : vec(ct_vec.vec) { }

// Assignment (copy-and-swap idiom)
void FandV_ct_vec::swap(FandV_ct_vec& a, FandV_ct_vec& b) {
  std::swap(a.vec, b.vec);
}
FandV_ct_vec& FandV_ct_vec::operator=(FandV_ct_vec ct_vec) {
  swap(*this, ct_vec);
  return(*this);
}

// Destructor
FandV_ct_vec::~FandV_ct_vec() {
  //Rcout << "DESTRUCT" << std::endl;
}

// Manipulate vector
void FandV_ct_vec::push(const FandV_ct& ct) {
  vec.push_back(ct);
}
void FandV_ct_vec::pushvec(const FandV_ct_vec& ct_vec) {
  for(unsigned int i=0; i<ct_vec.vec.size(); i++)
    vec.push_back(ct_vec.vec[i]);
}
void FandV_ct_vec::set(int i, const FandV_ct& ct_vec) {
  vec[i] = ct_vec;
}

// Access vector
int FandV_ct_vec::size() const {
  return(vec.size());
}
FandV_ct FandV_ct_vec::get(int i) const {
  return(vec[i]);
}
FandV_ct_vec FandV_ct_vec::subset(NumericVector i) const {
  FandV_ct_vec res;
  for(NumericVector::iterator it = i.begin(); it != i.end(); ++it) {
    res.push(vec[*it]);
  }
  return(res);
}
FandV_ct_vec FandV_ct_vec::without(NumericVector i) const { // NB must be sorted largest to smallest
  FandV_ct_vec res(vec);
  for(NumericVector::iterator it = i.begin(); it != i.end(); ++it) {
    res.vec.erase(res.vec.begin()+*it);
  }
  return(res);
}

// R level ops
FandV_ct_vec FandV_ct_vec::add(const FandV_ct_vec& x) const {
  int sz = vec.size(), xsz = x.vec.size();
  
  FandV_ct_vec res;
  if(sz>=xsz) {
    res.vec = vec;
    for(int i=0; i<sz; i++) {
      res.vec[i].addEq(x.vec[i%xsz]);
    }
  } else {
    res.vec = x.vec;
    for(int i=0; i<xsz; i++) {
      res.vec[i].addEq(vec[i%sz]);
    }
  }
  return(res);
}
FandV_ct_vec FandV_ct_vec::sub(const FandV_ct_vec& x) const {
  int sz = vec.size(), xsz = x.vec.size();
  
  FandV_ct_vec res;
  if(sz>=xsz) {
    res.vec = vec;
    for(int i=0; i<sz; i++) {
      res.vec[i] = vec[i].sub(x.vec[i%xsz]);
    }
  } else {
    res.vec = x.vec;
    for(int i=0; i<xsz; i++) {
      res.vec[i] = vec[i%sz].sub(x.vec[i]);
    }
  }
  return(res);
}
struct FandV_Mul : public Worker {   
  // Source vectors
  const std::vector<FandV_ct>* ctvec;
  const std::vector<FandV_ct>* ctvecX;
  const int sz, xsz;
  
  // Destination vector
  std::vector<FandV_ct>* res;
  
  // Constructors
  FandV_Mul(std::vector<FandV_ct>* res_, const std::vector<FandV_ct>* ctvec_, const std::vector<FandV_ct>* ctvecX_, int sz_, int xsz_) : sz(sz_), xsz(xsz_) { res = res_; ctvec = ctvec_; ctvecX = ctvecX_; }
  
  // Element wise multiply
  void operator()(std::size_t begin, std::size_t end) {
    for(; begin<end; begin++) {
      if(sz>=xsz)
        res->at(begin) = ctvec->at(begin).mul(ctvecX->at(begin%xsz));
      else
        res->at(begin) = ctvecX->at(begin).mul(ctvec->at(begin%sz));
    }
  }
};
FandV_ct_vec FandV_ct_vec::mulParallel(const FandV_ct_vec& x) const {
  int sz = vec.size(), xsz = x.vec.size();
  
  FandV_ct_vec res;
  if(sz>=xsz) {
    res.vec = vec;
  } else {
    res.vec = x.vec;
  }
  FandV_Mul mul(&(res.vec), &vec, &(x.vec), sz, xsz);
  parallelFor(0, res.vec.size(), mul);
  return(res);
}
FandV_ct_vec FandV_ct_vec::mulSerial(const FandV_ct_vec& x) const {
  int sz = vec.size(), xsz = x.vec.size();
  
  FandV_ct_vec res;
  if(sz>=xsz) {
    res.vec = vec;
    for(int i=0; i<sz; i++) {
      res.vec[i] = vec[i].mul(x.vec[i%xsz]);
    }
  } else {
    res.vec = x.vec;
    for(int i=0; i<xsz; i++) {
      res.vec[i] = x.vec[i].mul(vec[i%sz]);
    }
  }
  return(res);
}
FandV_ct_vec FandV_ct_vec::addct(const FandV_ct& ct) const {
  FandV_ct_vec res(vec);
  for(unsigned int i=0; i<vec.size(); i++) {
    res.vec[i].addEq(ct);
  }
  return(res);
}
FandV_ct_vec FandV_ct_vec::subct(const FandV_ct& ct, const int rev) const {
  FandV_ct_vec res(vec);
  for(unsigned int i=0; i<vec.size(); i++) {
    if(rev==0) {
      res.vec[i] = vec[i].sub(ct);
    } else {
      res.vec[i] = ct.sub(vec[i]);
    }
  }
  return(res);
}
struct FandV_MulCT : public Worker {   
  // Source vector
  const std::vector<FandV_ct>* ctvec;
  
  // Constant multiplier
  const FandV_ct* ct;
  
  // Destination vector
  std::vector<FandV_ct>* res;
  
  // Constructors
  FandV_MulCT(std::vector<FandV_ct>* res_, const std::vector<FandV_ct>* ctvec_, const FandV_ct* ct_)  { res = res_; ctvec = ctvec_; ct = ct_; }
  
  // Accumulate
  void operator()(std::size_t begin, std::size_t end) {
    for(; begin<end; begin++) {
      res->at(begin) = ctvec->at(begin).mul(*ct);
    }
  }
};
FandV_ct_vec FandV_ct_vec::mulctParallel(const FandV_ct& ct) const {
  FandV_ct_vec res(vec);
  FandV_MulCT mulct(&(res.vec), &vec, &ct);
  parallelFor(0, vec.size(), mulct);
  return(res);
}
FandV_ct_vec FandV_ct_vec::mulctSerial(const FandV_ct& ct) const {
  FandV_ct_vec res(vec);
  for(unsigned int i=0; i<vec.size(); i++) {
    res.vec[i] = vec[i].mul(ct);
  }
  return(res);
}

struct FandV_Sum : public Worker {   
  // Source vector
  const std::vector<FandV_ct>* input;
  
  // Accumulated value
  FandV_ct value;
  
  // Constructors
  FandV_Sum(const std::vector<FandV_ct>* input_) : value(input_->at(0).sub(input_->at(0))) { input = input_; }
  FandV_Sum(const FandV_Sum& sum, Split) : value(sum.input->at(0).sub(sum.input->at(0))) { input = sum.input; }
  
  // Accumulate
  void operator()(std::size_t begin, std::size_t end) {
    for(; begin<end; begin++) {
      value.addEq(input->at(begin));
    }
  }
  
  void join(const FandV_Sum& rhs) {
    value.addEq(rhs.value);
  }
};
FandV_ct FandV_ct_vec::sumParallel() const {
  FandV_Sum sum(&vec);
  parallelReduce(0, vec.size(), sum);
  return(sum.value);
}
FandV_ct FandV_ct_vec::sumSerial() const {
  FandV_ct res(vec[0]);
  
  for(unsigned int i=1; i<vec.size(); i++) {
    res.addEq(vec[i]);
  }
  
  return(res);
}

struct FandV_Prod : public Worker {   
  // Source vector
  const std::vector<FandV_ct>* input;
  
  // Accumulated value
  bool valueSet;
  FandV_ct value;
  
  // Constructors
  FandV_Prod(const std::vector<FandV_ct>* input_) : valueSet(false), value(input_->at(0).p, input_->at(0).rlkl, input_->at(0).rlki) { input = input_; }
  FandV_Prod(const FandV_Prod& prod, Split) : valueSet(false), value(prod.input->at(0).p, prod.input->at(0).rlkl, prod.input->at(0).rlki) { input = prod.input; }
  
  // Accumulate
  void operator()(std::size_t begin, std::size_t end) {
    for(; begin<end; begin++) {
      if(!valueSet) {
        value = input->at(begin);
        valueSet = true;
      } else {
        value = value.mul(input->at(begin));
      }
    }
  }
  
  void join(const FandV_Prod& rhs) {
    value = value.mul(rhs.value);
  }
};
FandV_ct FandV_ct_vec::prodParallel() const {
  FandV_Prod prod(&vec);
  parallelReduce(0, vec.size(), prod);
  return(prod.value);
}
FandV_ct FandV_ct_vec::prodSerial() const {
  FandV_ct res(vec[0]);
  
  for(unsigned int i=1; i<vec.size(); i++) {
    res = res.mul(vec[i]);
  }
  
  return(res);
}

struct FandV_InnerProd : public Worker {   
  // Source vector
  const std::vector<FandV_ct>* x;
  const std::vector<FandV_ct>* y;
  
  // Accumulated value
  bool resSet;
  FandV_ct res;
  
  // Constructors
  FandV_InnerProd(const std::vector<FandV_ct>* x_, const std::vector<FandV_ct>* y_) : resSet(false), res(x_->at(0).p, x_->at(0).rlkl, x_->at(0).rlki) { x = x_; y = y_; }
  FandV_InnerProd(const FandV_InnerProd& innerprod, Split) : resSet(false), res(innerprod.x->at(0).p, innerprod.x->at(0).rlkl, innerprod.x->at(0).rlki) { x = innerprod.x; y = innerprod.y; }
  
  // Accumulate
  void operator()(std::size_t begin, std::size_t end) {
    for(; begin<end; begin++) {
      if(!resSet) {
        res = x->at(begin).mul(y->at(begin));
        resSet = true;
      } else {
        res.addEq(x->at(begin).mul(y->at(begin)));
      }
    }
  }
  
  void join(const FandV_InnerProd& rhs) {
    res.addEq(rhs.res);
  }
};
FandV_ct FandV_ct_vec::innerprod(const FandV_ct_vec& x) const {
  FandV_InnerProd innerprod(&vec, &(x.vec));
  parallelReduce(0, vec.size(), innerprod);
  return(innerprod.res);
}

void FandV_ct_vec::show() const {
  Rcout << "Vector of " << vec.size() << " Fan and Vercauteren cipher texts\n";
}

// Save/load
void FandV_ct_vec::save(FILE* fp) const {
  fprintf(fp, "=> FHE package object <=\nRcpp_FandV_ct_vec\nn=%d\n", (int) vec.size());
  for(unsigned int i=0; i<vec.size(); i++) {
    vec[i].save(fp);
  }
}
FandV_ct_vec::FandV_ct_vec(FILE* fp, const FandV_par& p, FandV_rlk_locker* rlkl, size_t rlki) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE package object <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (CT_VEC)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_ct_vec\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a vector of ciphertext objects\n";
    free(buf);
    return;
  }
  
  int vecsz;
  len = fscanf(fp, "n=%d\n", &vecsz);
  for(int i=0; i<vecsz; i++) {
    FandV_ct ct(fp, p, rlkl, rlki);
    vec.push_back(ct);
  }
  free(buf);
}
