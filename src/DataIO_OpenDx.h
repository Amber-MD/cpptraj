#ifndef INC_DATAIO_OPENDX_H
#define INC_DATAIO_OPENDX_H
#include "DataIO.h"
/// Write OpenDx format data files.
class DataIO_OpenDx : public DataIO {
  public:
    DataIO_OpenDx() {}
    static DataIO* Alloc() { return (DataIO*)new DataIO_OpenDx(); }
    int ReadData(std::string const&, DataSetList&);
    int processWriteArgs(ArgList&)                 { return 0; }
    int WriteData(std::string const&,DataSetList const&,DimArray const&)         { return 1; }
    int WriteDataInverted(std::string const&,DataSetList const&,DimArray const&) { return 1; }
    int WriteData2D(std::string const&, DataSet const&, DimArray const&)         { return 1; }
    int WriteData3D(std::string const&, DataSet const&, DimArray const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    int LoadGrid(const char*, DataSet&);
};
#endif
