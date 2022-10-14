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
    /// Hold a list of pointers to DataSets in a DataSetList
    class SetPool;
    /// Hold a pointer to DataSet in a DataSetList and its original index
    class Set;
    typedef std::vector<Set> SetArray;
    /// Hold info for a NetCDF variable
    class NcVar;
    typedef std::vector<NcVar> VarArray;

    static void MarkVarIdRead(VarArray&, int);

    int read_cpptraj_vars(DataSetList&, std::string const&, VarArray&,
                          std::vector<unsigned int> const&);

    //int read_1d_var(DataSetList&, std::string const&, unsigned int, VarArray const&) const;

    int defineDim(std::string&, std::string const&,
                  unsigned int, std::string const&);

    int writeData_1D_xy(DataSet const*);

    int writeData_1D(DataSet const*, Dimension const&, SetArray const&);

    int writeData_2D(DataSet const*);

    int ncid_;                ///< Current netcdf ID
    int dimIdx_;              ///< Keep track of indices currently defined
    std::vector<int> varIDs_; ///< All variable IDs currently defined
    int* varIDptr_;           ///< Pointer to start of varIDs_
};
#endif
