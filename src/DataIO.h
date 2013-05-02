#ifndef INC_DATAIO_H
#define INC_DATAIO_H
#include "ArgList.h"
#include "DataSetList.h"
#include "CpptrajFile.h"
#include "Dimension.h"
// Class: DataIO
/// Base class that all DataIO objects inherit from.
class DataIO {
  public:
    virtual ~DataIO() {}
    typedef DataIO* (*AllocatorType)();
    /// Type to hold coordinate info for each dimension in DataSet.
    typedef std::vector<Dimension> DimArray;
    // ----- Inherited Functions -----------------
    virtual int ReadData(std::string const&, DataSetList&) = 0;
    virtual int processWriteArgs(ArgList&) = 0;
    virtual int WriteData(std::string const&, DataSetList const&, DimArray const&) = 0;
    virtual int WriteData2D(std::string const&, DataSet const&, DimArray const&) = 0;
    virtual int WriteData3D(std::string const&, DataSet const&, DimArray const&) = 0;
    virtual int WriteDataInverted(std::string const&, DataSetList const &, DimArray const&) = 0;
    virtual bool ID_DataFormat(CpptrajFile&);
  protected:
    // TODO: Move this to DataSet?
    static std::string SetupCoordFormat(size_t, Dimension const&, int, int);
};
#endif
