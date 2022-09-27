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
    typedef std::vector<Set> SetArray;

    int writeData_1D(DataSet const*, Dimension const&, SetArray const&);

    int ncid_;                ///< Current netcdf ID
    int dimIdx_;              ///< Keep track of indices currently defined
    std::vector<int> varIDs_; ///< All variable IDs currently defined
    int* varIDptr_;           ///< Pointer to start of varIDs_
};
#endif
