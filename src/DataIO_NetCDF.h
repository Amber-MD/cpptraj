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
    /// Hold info for a NetCDF dimension
    class NcDim;
    typedef std::vector<NcDim> DimArray;

    int readData_1D_xy(DataSet*, NcVar const&, VarArray&) const;

    int readData_2D(DataSet*, NcVar const&, VarArray&) const;

    int read_cpptraj_vars(DataSetList&, std::string const&, VarArray&) const;

    //int read_1d_var(DataSetList&, std::string const&, unsigned int, VarArray const&) const;

    int defineDim(std::string const&, unsigned int, std::string const&);

    NcVar defineVar(int, int, std::string const&, std::string const&, int) const;

    NcVar defineVar(int, int, std::string const&, std::string const&) const;

    int writeData_1D_xy(DataSet const*);

    int writeData_1D(DataSet const*, Dimension const&, SetArray const&);

    int writeData_2D(DataSet const*);

    int writeData_modes(DataSet const*);

    int ncid_;                 ///< Current netcdf ID
    bool user_specified_name_; ///< True if user specified a data set name on read
    DimArray Dimensions_;      ///< Array of all currently defined dimensions
};
#endif
