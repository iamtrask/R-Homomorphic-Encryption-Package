/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#ifndef FandV_ct_H
#define FandV_ct_H

#include "FandV_par.h"
#include "FandV_keys.h"

#include <fmpz_polyxx.h>
using namespace flint;

class FandV_ct {
  public:
    // Constructors
    FandV_ct(const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_);
    //FandV_ct(const FandV_par& p_, const FandV_rlk& rlk_);
    FandV_ct(const FandV_ct& ct);
    
    // Operators
    FandV_ct& operator=(FandV_ct ct);
    void swap(FandV_ct& a, FandV_ct& b);
    
    // R level ops
    FandV_ct add(const FandV_ct& c) const;
    void addEq(const FandV_ct& c); // += ... overwrites ct in place
    FandV_ct sub(const FandV_ct& c) const;
    FandV_ct mul(const FandV_ct& c) const;
    
    // Print out
    void show() const;
    
    // Save/load
    void save(FILE* fp) const;
    FandV_ct(FILE* fp, const FandV_par& p_, FandV_rlk_locker* rlkl_, size_t rlki_);
    
    // For performance keep public
    fmpz_polyxx c0, c1; // Polynomials
    FandV_par p;
    FandV_rlk_locker* rlkl;
    size_t rlki;
    int depth;
};

#endif
