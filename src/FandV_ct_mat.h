/*
 Louis Aslett (aslett@stats.ox.ac.uk)
 January 2015
*/

#ifndef FandV_ct_mat_H
#define FandV_ct_mat_H

#include <vector>

class FandV_ct;
class FandV_ct_vec;

class FandV_ct_mat {
  public:
    // Constructors
    FandV_ct_mat();
    FandV_ct_mat(const std::vector<FandV_ct>& v, const int nrow_, const int ncol_);
    FandV_ct_mat(const FandV_ct_mat& ct_mat);
    ~FandV_ct_mat();
    
    // Operators
    FandV_ct_mat& operator=(FandV_ct_mat ct_mat);
    void swap(FandV_ct_mat& a, FandV_ct_mat& b);
    
    // Manipulate vector
    //void push(const FandV_ct& ct);
    //void pushvec(const FandV_ct_vec& ct_vec);
    void set(int i, int j, const FandV_ct& ct);
    void setelt(int i, const FandV_ct& ct);
    void setmatrix(const FandV_ct_vec& ct_vec, int nrow_, int ncol_, int byrow);
    void reset(const FandV_ct& ct, const int nrow_, const int ncol_);
    
    // Access ...
    int size() const;
    // ... usual matrix
    FandV_ct get(int i) const; // Specify vector-like the element counting columnwise
    FandV_ct_mat subset(IntegerVector i, int nrow, int ncol) const; // vector indicies i chosen to form new matrix of nrow x ncol
    FandV_ct_vec subsetV(IntegerVector i) const;
    FandV_ct_mat t() const;
    
    // R level ops
    FandV_ct_mat add(const FandV_ct_mat& x) const;
    FandV_ct_mat mul(const FandV_ct_mat& x) const;
    FandV_ct_mat addct(const FandV_ct& ct) const;
    FandV_ct_mat mulctParallel(const FandV_ct& ct) const;
    FandV_ct_mat mulctSerial(const FandV_ct& ct) const;
    FandV_ct_mat mulctvecParallel(const FandV_ct_vec& ctvec) const;
    FandV_ct_mat mulctvecSerial(const FandV_ct_vec& ctvec) const;
    FandV_ct_mat matmulParallel(const FandV_ct_mat& y) const;
    FandV_ct_mat matmulSerial(const FandV_ct_mat& y) const;
    FandV_ct_mat TmatmulParallel(const FandV_ct_mat& y) const; // t(this) %*% y
    FandV_ct_mat matmulTParallel(const FandV_ct_mat& y) const; // this %*% t(y)
    FandV_ct_vec rowSumsParallel() const;
    FandV_ct_vec rowSumsSerial() const;
    FandV_ct_vec colSumsParallel() const;
    FandV_ct_vec colSumsSerial() const;
    // ADD VECTOR?  R DOES AND ADDS COLUMN WISE
    //FandV_ct sumParallel() const;
    //FandV_ct sumSerial() const;
    //FandV_ct prodParallel() const;
    //FandV_ct prodSerial() const;
    //FandV_ct innerprod(const FandV_ct_vec& x) const;
    
    // Print out
    void show() const;
    
    // Save/load
    void save(FILE* fp) const;
    FandV_ct_mat(FILE* fp, const FandV_par& p, FandV_rlk_locker* rlkl, size_t rlki);
    
    // For performance keep public
    int nrow;
    int ncol;
    std::vector<FandV_ct> mat;
};

#endif
