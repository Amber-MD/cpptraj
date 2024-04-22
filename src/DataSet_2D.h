#ifndef INC_DATASET_2D_H
#define INC_DATASET_2D_H
#include "DataSet.h"
/// Interface for 2D DataSets.
class DataSet_2D : public DataSet {
  public:
    /// Kind of matrix.
    enum MatrixKindType { FULL = 0, HALF, TRI };
    DataSet_2D() {}
    DataSet_2D(DataSet::DataType tIn, TextFormat const& fIn) : 
      DataSet(tIn, MATRIX_2D, fIn, 2) {}
    // TODO enable Append?
    int Append(DataSet*) { return 1; }
    // ----- DataSet functions -------------------
    int Allocate(SizeArray const& s) { return Allocate2D(s[0], s[1]); }
    // ----- DataSet_2D functions ----------------
    // TODO: 1 allocate using MatrixKind?
    /// Set up matrix for given # columns and rows. // TODO replace with Allocate
    virtual int Allocate2D(size_t, size_t) = 0;
    /// Set up symmetric matrix with diagonal.
    virtual int AllocateHalf(size_t) = 0;
    /// Set up symmetrix matrix with no diagonal.
    virtual int AllocateTriangle(size_t) = 0;

    /// Add given value to element at specified col/row
    virtual void UpdateElement(size_t, size_t, double) = 0;
    /// Add given value to element at specified 1D index
    virtual void UpdateElement(size_t, double) = 0;
    /// Set element at specified col/row to given value
    virtual void SetElement(size_t, size_t, double) = 0;
    /// Set element at specified 1D index to given value
    virtual void SetElement(size_t, double) = 0;

    /// \return Data from matrix at col/row 
    virtual double GetElement(size_t, size_t) const = 0;
    /// \return Data from underlying matrix array.
    virtual double GetElement(size_t) const = 0;
    /// \return the number of rows.
    virtual size_t Nrows() const = 0;
    /// \return the number of columns.
    virtual size_t Ncols() const = 0;
    /// \return double array containing matrix elements.
    virtual double* MatrixArray() const = 0;
    /// \return the kind of matrix, full/half/triangle.
    virtual MatrixKindType MatrixKind() const = 0;
    /// \return const pointer to underlying matrix.
    virtual const void* MatrixPtr() const = 0;
    /// \return pointer to underlying matrix.
    virtual void* MatrixPtr() = 0;
    /// Clear matrix
    virtual void Clear() = 0;
    /// Normalize matrix
    virtual void Normalize(double) = 0;
    // -------------------------------------------
    enum RetType { OK = 0, DIM_MISMATCH, ALLOC_ERR, ERR };
    /// Set this matrix with the result of M1 x M2
    RetType Multiply(DataSet_2D const& M1, DataSet_2D const& M2) {
      if (M1.Ncols() != M2.Nrows()) return DIM_MISMATCH;
      unsigned int len = M1.Ncols();
      if (Allocate2D( M2.Ncols(), M1.Nrows() )) return ALLOC_ERR;
      unsigned int idx = 0;
      for (unsigned int row = 0; row != Nrows(); row++) {
        for (unsigned int col = 0; col != Ncols(); col++) {
          double sum = 0;
          for (unsigned int k = 0; k != len; k++)
            sum += M1.GetElement(k, row) * M2.GetElement(col, k);
          SetElement(idx++, sum);
        }
      }
      return OK;
    }
    /// Set this matrix with the result of M1 x M2^T
    RetType Multiply_M2transpose(DataSet_2D const& M1, DataSet_2D const& M2) {
      if (M1.Ncols() != M2.Ncols()) return DIM_MISMATCH;
      unsigned int len = M1.Ncols();
      if (Allocate2D( M2.Nrows(), M1.Nrows() )) return ALLOC_ERR;
      unsigned int idx = 0;
      for (unsigned int row = 0; row != Nrows(); row++) {
        for (unsigned int col = 0; col != Ncols(); col++) {
          double sum = 0;
          for (unsigned int k = 0; k != len; k++)
            sum += M1.GetElement(k, row) * M2.GetElement(k, col);
          SetElement(idx++, sum);
        }
      }
      return OK;
    }
    /// Set this matrix with the result of M1^T x M2
    RetType Multiply_M1transpose(DataSet_2D const& M1, DataSet_2D const& M2) {
      if (M1.Nrows() != M2.Nrows()) return DIM_MISMATCH;
      unsigned int len = M1.Nrows();
      if (Allocate2D( M2.Ncols(), M1.Ncols() )) return ALLOC_ERR;
      unsigned int idx = 0;
      for (unsigned int row = 0; row != Nrows(); row++) {
        for (unsigned int col = 0; col != Ncols(); col++) {
          double sum = 0;
          for (unsigned int k = 0; k != len; k++)
            sum += M1.GetElement(row, k) * M2.GetElement(col, k);
          SetElement(idx++, sum);
        }
      }
      return OK;
    }

    /// Set this matrix with the result of M1^T
    RetType TransposeOf(DataSet_2D const& M1) {
      RetType ret = ERR;
      MatrixKindType m1kind = M1.MatrixKind();
      if (m1kind == FULL) {
        if (Allocate2D( M1.Nrows(), M1.Ncols() )) return ALLOC_ERR;
        ret = OK;
        for (unsigned int row = 0; row < M1.Nrows(); row++) {
          for (unsigned int col = 0; col < M1.Ncols(); col++)
            SetElement( row, col, M1.GetElement(col, row) );
        }
      } else if (m1kind == HALF) {
        if (AllocateHalf( M1.Ncols() )) return ALLOC_ERR;
        ret = OK;
        for (unsigned int row = 0; row < M1.Nrows(); row++) {
          for (unsigned int col = row; col < M1.Ncols(); col++)
            SetElement( row, col, M1.GetElement(col, row) );
        }
      } else if (m1kind == TRI) {
        if (AllocateTriangle( M1.Ncols() )) return ALLOC_ERR;
        ret = OK;
        for (unsigned int row = 0; row < M1.Nrows(); row++) {
          for (unsigned int col = row + 1; col < M1.Ncols(); col++)
            SetElement( row, col, M1.GetElement(col, row) );
        }
      }
      return ret;
    }
    /// \return True if matrix is symmetric
    bool IsSymmetric() const;
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
};
#endif
