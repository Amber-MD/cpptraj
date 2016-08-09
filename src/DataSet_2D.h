#ifndef INC_DATASET_2D_H
#define INC_DATASET_2D_H
#include "DataSet.h"
#include "CpptrajFile.h"
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
    int Allocate(SizeArray const& s) { return Allocate2D(s[0], s[1]); } // FIXME bounds check
    // TODO: 1 allocate using MatrixKind?
    /// Set up matrix for given # columns and rows. // TODO replace with Allocate
    virtual int Allocate2D(size_t, size_t) = 0;
    /// Set up symmetric matrix with diagonal.
    virtual int AllocateHalf(size_t) = 0;
    /// Set up symmetrix matrix with no diagonal.
    virtual int AllocateTriangle(size_t) = 0;
    /// Add given value to element at specified col/row
    virtual void UpdateElement(size_t, size_t, double) = 0;
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
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
};
#endif
