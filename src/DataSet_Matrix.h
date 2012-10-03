#ifndef INC_DATASET_MATRIX_H
#define INC_DATASET_MATRIX_H
#include "DataSet.h"
/// Matrix, hold average over frames
class DataSet_Matrix : public DataSet {
  public:
    enum Matrix_Type {
      MATRIX_NULL=0, MATRIX_DIST,      MATRIX_COVAR, MATRIX_MWCOVAR,
      MATRIX_CORREL, MATRIX_DISTCOVAR, MATRIX_IDEA,  MATRIX_IRED
    };

    DataSet_Matrix();
    ~DataSet_Matrix();

    void SetType( Matrix_Type typeIn ) { type_ = typeIn; }

    // Iterator over matrix elements
    class iterator : public std::iterator<std::forward_iterator_tag, double>
    {
      public:
        iterator() : ptr_(0) {}
        iterator(const iterator& rhs) : ptr_(rhs.ptr_) {}
        iterator(double* pin) : ptr_(pin) {}
        iterator& operator=(const iterator& rhs) {
          if (this == &rhs) return *this;
          ptr_ = rhs.ptr_;
          return *this;
        }
        // Relations
        bool operator==(const iterator& rhs) { return (ptr_==rhs.ptr_);}
        bool operator!=(const iterator& rhs) { return (ptr_!=rhs.ptr_);}
        // Increment
        iterator& operator++() {
          ++ptr_;
          return *this;
        }
        iterator operator++(int) {
          iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        // Value
        double& operator*() { return *ptr_; }
        // Address
        double* operator->() { return ptr_; }
      private:
        double* ptr_;
    };
    iterator begin() { return mat_; }
    iterator end() { return (mat_ + matsize_); }
  private:
    double* mat_;     ///< Hold matrix elements
    double* vect_;    ///< Hold diagonal elements | avg coords
    double* vect2_;   ///< Square of vect_. May not need to be stored.
    double* mass_;    ///< Hold masses. Currently only for square (i.e. size is nrows)?
    int matsize_;     ///< Total number of matrix elements
    int nrows_;       ///< Number of rows in the matrix
    int ncols_;       ///< Number of columns in the matrix
    int vectsize_;    ///< Sizes of vect | vect2
    int snap_;        ///< Number of snapshots
    Matrix_Type type_; ///< Type of matrix.
};
#endif
