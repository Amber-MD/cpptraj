#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
// Class: DataIO_Std
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();
    static DataIO* Alloc() { return (DataIO*)new DataIO_Std(); }
    int ReadData(std::string const&,DataSetList&);
    int processWriteArgs(ArgList&);
    int WriteData(std::string const&,DataSetList const&,DimArray const&);
    int WriteDataInverted(std::string const&,DataSetList const&,DimArray const&);
    int WriteData2D(std::string const&, DataSet const&, DimArray const&);
    int WriteData3D(std::string const&, DataSet const&, DimArray const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    bool hasXcolumn_;
    bool writeHeader_;
    bool square2d_;
};
#endif
