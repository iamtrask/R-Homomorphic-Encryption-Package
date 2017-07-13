/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 August 2014
*/

#ifndef FandV_ct_vec_H
#define FandV_ct_vec_H

class FandV_ct;

class FandV_ct_vec {
  public:
    // Constructors
    FandV_ct_vec();
    FandV_ct_vec(const std::vector<FandV_ct> v);
    FandV_ct_vec(const FandV_ct_vec& ct_vec);
    ~FandV_ct_vec();
    
    // Operators
    FandV_ct_vec& operator=(FandV_ct_vec ct_vec);
    void swap(FandV_ct_vec& a, FandV_ct_vec& b);
    
    // Manipulate vector
    void push(const FandV_ct& ct);
    void pushvec(const FandV_ct_vec& ct_vec);
    void set(int i, const FandV_ct& ct_vec);
    
    // Access vector
    int size() const;
    FandV_ct get(int i) const;
    FandV_ct_vec subset(NumericVector i) const;
    FandV_ct_vec without(NumericVector i) const; // NB must be sorted largest to smallest
    
    // R level ops
    FandV_ct_vec add(const FandV_ct_vec& x) const;
    FandV_ct_vec sub(const FandV_ct_vec& x) const;
    FandV_ct_vec mulParallel(const FandV_ct_vec& x) const;
    FandV_ct_vec mulSerial(const FandV_ct_vec& x) const;
    FandV_ct_vec addct(const FandV_ct& ct) const;
    FandV_ct_vec subct(const FandV_ct& ct, const int rev) const;
    FandV_ct_vec mulctParallel(const FandV_ct& ct) const;
    FandV_ct_vec mulctSerial(const FandV_ct& ct) const;
    FandV_ct sumParallel() const;
    FandV_ct sumSerial() const;
    FandV_ct prodParallel() const;
    FandV_ct prodSerial() const;
    FandV_ct innerprod(const FandV_ct_vec& x) const;
    
    // Print out
    void show() const;
    
    // Save/load
    void save(FILE* fp) const;
    FandV_ct_vec(FILE* fp, const FandV_par& p, FandV_rlk_locker* rlkl, size_t rlki);
    
    // For performance keep public
    std::vector<FandV_ct> vec;
};

#endif
