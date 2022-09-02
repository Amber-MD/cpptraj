#ifndef INC_DATAIO_NETCDF_H
#define INC_DATAIO_NETCDF_H
#include "DataIO.h"
/// Generic NetCDF DataSet 
class DataIO_NetCDF : public DataIO {
  public:
    DataIO_NetCDF();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_NetCDF(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    class SetPool;
    class Set;
};
#endif
