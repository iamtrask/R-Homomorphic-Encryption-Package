/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#ifndef FandV_keys_H
#define FandV_keys_H

#include <Rcpp.h>
using namespace Rcpp;

#include <fmpz_polyxx.h>
using namespace flint;

#include "FandV_par.h"

#include <vector>

class FandV_ct;
class FandV_ct_vec;
class FandV_ct_mat;
class FandV_sk;
class FandV_pk;

class FandV_rlk {
  public:
    // Constructors
    FandV_rlk();
    FandV_rlk(const FandV_rlk& rlk);
    
    // Relinearise
    //int relin(FandV_ct& ct);
    
    // Print
    void show();
    
    friend void FandV_par::keygen(FandV_pk& pk, FandV_sk& sk, FandV_rlk& rlk);
    
    // Save/load
    void save(FILE* fp) const;
    FandV_rlk(FILE* fp);
    
    fmpz_polyxx rlk00, rlk01, rlk10, rlk11;
};

class FandV_rlk_locker {
  public:
    // Constructors
    FandV_rlk_locker();
    
    // Add a relin key to the locker and return the index
    int add(const FandV_rlk &rlk);
    void show() const;
    
    // The locker containing relin keys
    std::vector<FandV_rlk> x;
};

class FandV_pk {
  public:
    // Constructors/Destructors
    FandV_pk(FandV_rlk_locker* rlkl, size_t rlki);
    FandV_pk(const FandV_pk& pk);
    
    // Encrypt
    void enc(int m, FandV_ct& ct) const;
    void encvec(IntegerVector m, FandV_ct_vec& ctvec);
    void encmat(IntegerVector m, int nrow, int ncol, FandV_ct_mat& ctmat);
    
    // Print
    void show();
    
    friend void FandV_par::keygen(FandV_pk& pk, FandV_sk& sk, FandV_rlk& rlk);

    // Save/load
    void save(FILE* fp) const;
    FandV_pk(FILE* fp, const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_);

    FandV_par p;
    FandV_rlk_locker* rlkl;
    size_t rlki;
    
  private:
    fmpz_polyxx p0, p1; // Cyclotomic polynomial defining ring modulo
};

class FandV_sk {
  public:
    // Constructors
    FandV_sk();
    FandV_sk(const FandV_sk& sk);
    
    // Decrypt
    std::vector<int> decraw(const FandV_ct& ct) const;
    std::string dec(const FandV_ct& ct) const;
    
    // Print
    void show();
    
    friend void FandV_par::keygen(FandV_pk& pk, FandV_sk& sk, FandV_rlk& rlk);
    
    // Save/load
    void save(FILE* fp) const;
    FandV_sk(FILE* fp);

  private:
    fmpz_polyxx s; // Cyclotomic polynomial defining ring modulo
};

#endif
