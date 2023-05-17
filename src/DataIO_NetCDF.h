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
#ifdef BINTRAJ
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
    /// Hold array of integers
//    typedef std::vector<int> Iarray;

    std::vector<Dimension> getVarIndexInfo(int&, int, int) const;

    int get_1D_var_dimlen(size_t&, int, const char*) const;

    inline unsigned int dimLen(int) const;

    int readData_1D_xy(DataSet*, NcVar const&, VarArray&) const;

    int readData_1D_string(DataSet*, NcVar const&, VarArray&) const;

    int readData_1D_unsignedInt(DataSet*, NcVar const&, VarArray&) const;

    int readData_Mat3x3(DataSet*, NcVar const&, VarArray&) const;

    int readData_1D_vector(DataSet*, NcVar const&, VarArray&) const;

    int readData_1D_vector_scalar(DataSet*, NcVar const&, VarArray&) const;

    int readData_1D(DataSet*, NcVar const&, VarArray&) const;

    int readData_cluster_pwmatrix(DataSet*, NcVar const&, VarArray&) const;

    int readData_2D(DataSet*, NcVar const&, VarArray&) const;

    int readData_3D(DataSet*, NcVar const&, VarArray&) const;

    int readData_modes(DataSet*, NcVar const&, VarArray&) const;

    int read_cpptraj_vars(DataSetList&, std::string const&, VarArray&) const;

    //int read_1d_var(DataSetList&, std::string const&, unsigned int, VarArray const&) const;

    int defineDim(std::string const&, unsigned int, std::string const&);

    NcVar defineVar(std::vector<int> const&, int, std::string const&, std::string const&, int) const;

    NcVar defineVar(int, int, std::string const&, std::string const&, int) const;

    NcVar defineVar(int, int, std::string const&, std::string const&) const;

    int writeData_vector_scalar(DataSet const*);

    int writeData_1D(DataSet const*, Dimension const&, SetArray const&);

    int writeData_2D(DataSet const*);

    int writeData_3D(DataSet const*);

    int writeData_cluster_pwmatrix(DataSet const*);

    int writeData_modes(DataSet const*);

    int ncid_;                        ///< Current netcdf ID
    std::string user_specified_name_; ///< Set if user specified a data set name on read
    DimArray Dimensions_;             ///< Array of all currently defined dimensions
#   endif /* BINTRAJ */
};
#endif
