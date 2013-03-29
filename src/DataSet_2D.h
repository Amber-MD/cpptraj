#ifndef INC_DATASET_2D_H
#define INC_DATASET_2D_H
#include "DataSet.h"
#include "CpptrajFile.h"
/// Interface for 2D DataSets.
class DataSet_2D : public DataSet {
  public:
    enum MatrixType {
      NO_OP=0, DIST, COVAR, MWCOVAR, CORREL, DISTCOVAR, IDEA, IRED, NMAT
    };

    DataSet_2D() : type_(NO_OP) {}
    DataSet_2D(DataSet::DataType tIn, int wIn, int pIn) : 
      DataSet(tIn, wIn, pIn, 2), type_(NO_OP) {}
    /// Set up matrix for given # rows and columns.
    virtual int Allocate2D(size_t, size_t) = 0;
    virtual int AllocateHalf(size_t) = 0;
    virtual int AllocateTriangle(size_t) = 0;
    /// Write 2D data to file (2D)
    virtual void Write2D(CpptrajFile&,int,int) const = 0;
    /// \return Data from matrix at row/col
    virtual double GetElement(size_t, size_t) const = 0;
    /// \return the number of rows.
    virtual size_t Nrows() const = 0;
    /// \return the number of columns.
    virtual size_t Ncols() const = 0;
    /// \return double array containing matrix elements.
    virtual double* MatrixArray() const = 0;

    static const char* MatrixTypeString[];
    static const char* MatrixOutputString[];
  private:
    MatrixType type_;
};
#endif
