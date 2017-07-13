/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "FandV_ct.h"
#include "FandV.h"

// Construct from parameters
FandV_ct::FandV_ct(const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_) : p(p_), rlkl(rlkl_), rlki(rlki_), depth(0) { }

// Copy constructor
FandV_ct::FandV_ct(const FandV_ct& ct) : c0(ct.c0), c1(ct.c1), p(ct.p), rlkl(ct.rlkl), rlki(ct.rlki), depth(ct.depth) { }

// Assignment (copy-and-swap idiom)
void FandV_ct::swap(FandV_ct& a, FandV_ct& b) {
  std::swap(a.c0, b.c0);
  std::swap(a.c1, b.c1);
  std::swap(a.p, b.p);
  std::swap(a.rlkl, b.rlkl);
  std::swap(a.rlki, b.rlki);
  std::swap(a.depth, b.depth);
}
FandV_ct& FandV_ct::operator=(FandV_ct ct) {
  swap(*this, ct);
  return(*this);
}

// R level ops
FandV_ct FandV_ct::add(const FandV_ct& c) const {
  FandV_ct res(p, rlkl, rlki);
  res.depth = std::max(depth, c.depth);
  
  res.c0 = c0+c.c0;
  res.c1 = c1+c.c1;
  
  return(res);
}
void FandV_ct::addEq(const FandV_ct& c) {  
  depth = std::max(depth, c.depth);
  
  c0 += c.c0;
  c1 += c.c1;
}

FandV_ct FandV_ct::sub(const FandV_ct& c) const {
  FandV_ct res(p, rlkl, rlki);
  res.depth = std::max(depth, c.depth);
  
  res.c0 = c0-c.c0;
  res.c1 = c1-c.c1;
  
  return(res);
}

FandV_ct FandV_ct::mul(const FandV_ct& c) const {
  fmpz_polyxx c2, res2;
  c2.realloc(p.Phi.length());
  res2.realloc(p.Phi.length());
  fmpzxx one(1);
  FandV_ct res(p, rlkl, rlki);
  res.depth = depth+c.depth+1;
  
  // c0
  //res.c0 = ((c0*c.c0)%p.Phi); Rcout << res.c0 << "\n"; // Following indented lines are 2x faster at doing modulo cyclotomic poly
    res.c0 = c0*c.c0;
    for(int i=0; i<p.Phi.length()-1; i++) {
      res.c0.set_coeff(i, res.c0.get_coeff(i)-res.c0.get_coeff(i+p.Phi.length()-1));
      res.c0.set_coeff(i+p.Phi.length()-1, 0);
    }
    res.c0.set_coeff(2*p.Phi.length()-2, 0);
    
  res2   = (p.t*res.c0)%p.q;
  res.c0 = (p.t*res.c0)/p.q;
  for(int i=0; i<p.Phi.length(); i++) {
    if(res2.get_coeff(i) > p.q/2)
      res.c0.set_coeff(i, res.c0.get_coeff(i)+one);
  }
  fmpz_polyxx_q(res.c0, p.q);
  
  
  // c1
  //res.c1 = ((c0*c.c1 + c1*c.c0)%p.Phi); // Following indented lines are 2x faster at doing modulo cyclotomic poly
    res.c1 = c0*c.c1 + c1*c.c0;
    for(int i=0; i<p.Phi.length()-1; i++) {
      res.c1.set_coeff(i, res.c1.get_coeff(i)-res.c1.get_coeff(i+p.Phi.length()-1));
      res.c1.set_coeff(i+p.Phi.length()-1, 0);
    }
    res.c1.set_coeff(2*p.Phi.length()-2, 0);
  
  res2   = (p.t*res.c1)%p.q;
  res.c1 = (p.t*res.c1)/p.q;
  for(int i=0; i<p.Phi.length(); i++) {
    if(res2.get_coeff(i) > p.q/2)
      res.c1.set_coeff(i, res.c1.get_coeff(i)+one);
  }
  fmpz_polyxx_q(res.c1, p.q);
  
  
  // c2
  //c2 = ((c1*c.c1)%p.Phi); // Following indented lines are 2x faster at doing modulo cyclotomic poly
    c2 = c1*c.c1;
    for(int i=0; i<p.Phi.length()-1; i++) {
      c2.set_coeff(i, c2.get_coeff(i)-c2.get_coeff(i+p.Phi.length()-1));
      c2.set_coeff(i+p.Phi.length()-1, 0);
    }
    c2.set_coeff(2*p.Phi.length()-2, 0);
  
  res2 = (p.t*c2)%p.q;
  c2 = (p.t*c2)/p.q;
  for(int i=0; i<p.Phi.length(); i++) {
    if(res2.get_coeff(i) > p.q/2)
      c2.set_coeff(i, c2.get_coeff(i)+one);
  }
  fmpz_polyxx_q(c2, p.q);
  
  
  // relin
  for(int i=0; i<p.Phi.length(); i++) {
    res2.set_coeff(i, c2.get_coeff(i)%p.T);
    c2.set_coeff(i, c2.get_coeff(i)/p.T);
  }
  
  FandV_rlk& rlk = (rlkl->x)[rlki];
  //res.c0 = res.c0 + ((rlk.rlk00*res2)%p.Phi) + ((rlk.rlk10*c2)%p.Phi); // Following indented lines are 2x faster at doing modulo cyclotomic poly
    res.c0 = res.c0 + rlk.rlk00*res2 + rlk.rlk10*c2;
    for(int i=0; i<p.Phi.length()-1; i++) {
      res.c0.set_coeff(i, res.c0.get_coeff(i)-res.c0.get_coeff(i+p.Phi.length()-1));
      res.c0.set_coeff(i+p.Phi.length()-1, 0);
    }
    res.c0.set_coeff(2*p.Phi.length()-2, 0);
    
  fmpz_polyxx_q(res.c0, p.q);
  
  //res.c1 = res.c1 + ((rlk.rlk01*res2)%p.Phi) + ((rlk.rlk11*c2)%p.Phi); // Following indented lines are 2x faster at doing modulo cyclotomic poly
    res.c1 = res.c1 + rlk.rlk01*res2 + rlk.rlk11*c2;
    for(int i=0; i<p.Phi.length()-1; i++) {
      res.c1.set_coeff(i, res.c1.get_coeff(i)-res.c1.get_coeff(i+p.Phi.length()-1));
      res.c1.set_coeff(i+p.Phi.length()-1, 0);
    }
    res.c1.set_coeff(2*p.Phi.length()-2, 0);
  
  fmpz_polyxx_q(res.c1, p.q);
  
  return(res);
}

void FandV_ct::show() const {
  Rcout << "Fan and Vercauteren cipher text\n";
  Rcout << "( c\u2080 = ";
  printPoly(c0);
  Rcout << ",\nc\u2081 = ";
  printPoly(c1);
  Rcout << " )\n";
}

// Save/load
void FandV_ct::save(FILE* fp) const {
  fprintf(fp, "=> FHE pkg obj <=\nRcpp_FandV_ct\n");
  // c0
  print(fp, c0);
  fprintf(fp, "\n");
  // c1
  print(fp, c1);
  fprintf(fp, "\n");
  // depth
  fprintf(fp, "%d\n", depth);
}
FandV_ct::FandV_ct(FILE* fp, const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_) : p(p_), rlkl(rlkl_), rlki(rlki_) {
  // Check for header line
  char *buf = NULL; size_t bufn = 0;
  size_t len;
  len = getline(&buf, &bufn, fp);
  if(strncmp("=> FHE pkg obj <=\n", buf, len) != 0) {
    Rcout << "Error: file does not contain an FHE object (CT)\n";
    free(buf);
    return;
  }
  len = getline(&buf, &bufn, fp);
  if(strncmp("Rcpp_FandV_ct\n", buf, len) != 0) {
    Rcout << "Error: file does not contain a single ciphertext object\n";
    free(buf);
    return;
  }
  
  // c0
  read(fp, c0);
  // c1
  read(fp, c1);
  // depth
  len = fscanf(fp, "%d\n", &depth);
  
  free(buf);
}
