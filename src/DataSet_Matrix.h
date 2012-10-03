#ifndef INC_DATASET_MATRIX_H
#define INC_DATASET_MATRIX_H
#include "DataSet.h"
/// Matrix, hold average over frames
class DataSet_Matrix : public DataSet {
  public:
    enum matrixMode {
      MATRIX_NULL=0, MATRIX_DIST,      MATRIX_COVAR, MATRIX_MWCOVAR,
      MATRIX_CORREL, MATRIX_DISTCOVAR, MATRIX_IDEA,  MATRIX_IRED
    };

    DataSet_Matrix();
    ~DataSet_Matrix();
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
    matrixMode type_; ///< Type of matrix.
};
#endif
