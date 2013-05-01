#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include "DataSet.h"
#include "CpptrajFile.h"
/// Interface for 3D DataSets.
class DataSet_3D : public DataSet {
  public:
    DataSet_3D() {}
    DataSet_3D(DataSet::DataType tIn, int wIn, int pIn) :
      DataSet(tIn, wIn, pIn, 3) {}
    /// Set up grid for given # x, y, and z points.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;
    /// Write 3D data to file.
    virtual void Write3D(CpptrajFile&,int,int,int) const = 0;
    /// \return Data from grid at x/y/z point.
    virtual double GetElement(size_t, size_t, size_t) const = 0;
    /// \return size of X dimension.
    virtual size_t NX() const = 0;
    /// \return size of Y dimension.
    virtual size_t NY() const = 0;
    /// \return size of Z dimension.
    virtual size_t NZ() const = 0;
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
};
#endif
