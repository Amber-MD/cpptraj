#ifndef INC_DATASET_1D_H
#define INC_DATASET_1D_H
#include "DataSet.h"
#include "CpptrajFile.h"
/// Class that all 1D DataSets will inherit.
class DataSet_1D : public DataSet {
  public:
    DataSet_1D() {}
    DataSet_1D(DataSet::DataType tIn, int wIn, int pIn) : DataSet(tIn, wIn, pIn, 1) {}
    virtual ~DataSet_1D() {}
    /// Allocate memory for a certain number of frames (1D).
    virtual int Allocate1D(size_t) = 0;
    /// Add data to the DataSet.
    /** A pointer to the data is passed in as void - it is up to the 
      * inheriting class to cast it. The X value for the data is passed 
      * in as well. It is expected that each successive X value will
      * be greater than the preceeding one (does not need to be
      * consecutive however). 
      */
    virtual void Add( size_t, const void* ) = 0;
    /// Write data at frame to file (1D)
    virtual void WriteBuffer(CpptrajFile&,size_t) const = 0;
    /// \return data from data set as double precision (1D)
    virtual double Dval(size_t) const = 0;
};
#endif
