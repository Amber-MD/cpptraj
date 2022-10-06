#include "DataIO_NetCDF.h"
#include "CpptrajStdio.h"
#ifdef BINTRAJ
# include <map>
# include <string>
# include <netcdf.h>
# include "DataSet_1D.h"
# include "NC_Routines.h"
# include "StringRoutines.h"
# include "Version.h"
#endif

/// CONSTRUCTOR
DataIO_NetCDF::DataIO_NetCDF() :
  DataIO(true, true, true), // Valid for 1D, 2D, 3D
  ncid_(-1),
  dimIdx_(0),
  varIDptr_(0)
{

}

// DataIO_NetCDF::ID_DataFormat()
bool DataIO_NetCDF::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  // Read first 3 bytes
  unsigned char magic[3];
  magic[0] = 0;
  magic[1] = 0;
  magic[2] = 0;
  infile.Read(magic, 3);
  infile.CloseFile();
  if (magic[0] == 'C' && magic[1] == 'D' && magic[2] == 'F') {
#   ifdef BINTRAJ
    return true;
#   else
    mprintf("Warning: '%s' is a NetCDF file but CPPTRAJ was compiled without NetCDF support.\n",
            infile.Filename().full());
#   endif
  }
  return false;
}

// DataIO_NetCDF::ReadHelp()
void DataIO_NetCDF::ReadHelp()
{

}

// DataIO_NetCDF::processReadArgs()
int DataIO_NetCDF::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_NetCDF::ReadData()
int DataIO_NetCDF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
# ifdef BINTRAJ
  if (NC::CheckErr( nc_open( fname.full(), NC_NOWRITE, &ncid_ ) != NC_NOERR )) {
    mprinterr("Error: Could not open NetCDF data file '%s'\n", fname.full());
    return 1;
  }

  // Get the number of dimensions, variables, attributes, and ID of the
  // unlimited dimension (if any).
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  if (NC::CheckErr( nc_inq(ncid_, &ndimsp, &nvarsp, &ngattsp, &unlimdimidp) )) {
    mprinterr("Error: Could not get NetCDF data file information.\n");
    return 1;
  }
  mprintf("DEBUG: '%s' : ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
          fname.full(), ndimsp, nvarsp, ngattsp, unlimdimidp);

  nc_close( ncid_ );
  return 0;
# else  
  return 1;
# endif
}

// -----------------------------------------------------------------------------
// DataIO_NetCDF::WriteHelp()
void DataIO_NetCDF::WriteHelp()
{

}

// DataIO_NetCDF::processWriteArgs()
int DataIO_NetCDF::processWriteArgs(ArgList& argIn)
{

  return 0;
}

/** Hold a pool of pointers to DataSets in the list. They will be marked off
  * as they are used.
  */
class DataIO_NetCDF::SetPool {
  public:
    /// CONSTRUCTOR - place sets from DataSetList in this pool
    SetPool(DataSetList const& dsl) : nUsed_(0) {
      sets_.reserve( dsl.size() );
      isUsed_.assign( dsl.size(), false );
      for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
        sets_.push_back( *it );
    }
    /// \return set at idx
    DataSet const* Set(unsigned int idx) const { return sets_[idx]; }
    /// \return Number of sets
    unsigned int Nsets() const { return sets_.size(); }
    /// \return true if set at idx has been used
    bool IsUsed(unsigned int idx) const { return isUsed_[idx]; }
    /// \return true if all sets have been marked as used
    bool AllUsed() const { return (nUsed_ == isUsed_.size()); }

    /// Mark set at idx as used
    void MarkUsed(unsigned int idx) {
      isUsed_[idx] = true;
      nUsed_++;
    }
  private:
    std::vector<DataSet const*> sets_;
    std::vector<bool> isUsed_;
    unsigned int nUsed_;
};

/** Hold a pointer to DataSet and its original index. Used to refer back to
  * original DataSetList from sets in a SetPool.
  */
class DataIO_NetCDF::Set {
  public:
    Set(DataSet const* ds, int oidx) : ds_(ds), oidx_(oidx) {}

    DataSet const* DS() const { return ds_; }
    int OriginalIdx() const { return oidx_; }
  private:
    DataSet const* ds_;
    int oidx_;
};

#ifdef BINTRAJ
/// Tell netcdf file to end define mode
static inline int EndDefineMode(int ncid) {
  // End netcdf definitions
  if (NC::CheckErr(nc_enddef(ncid))) {
    mprinterr("NetCDF data error on ending definitions.");
    return 1;
  }
  return 0;
}

/// Tell netcdf file to enter define mode
static inline int EnterDefineMode(int ncid) {
  int err = nc_redef(ncid);
  if (err != NC_NOERR && err != NC_EINDEFINE) {
    NC::CheckErr( err );
    return 1;
  }
  return 0;
}

/// Add a string as text attribute to a variable
static inline int AddDataSetStringAtt(std::string const& str, const char* desc, int ncid, int varid)
{
  if (!str.empty()) {
    if (NC::CheckErr(nc_put_att_text(ncid, varid, desc, str.size(), str.c_str()) ))
    {  
      mprinterr("Error: Writing variable %s.\n", desc);
      return 1;
    }
  }
  return 0;
}

/// Add an integer as integer attribute
static inline int AddDataSetIntAtt(int ival, const char* desc, int ncid, int varid)
{
  // NOTE: -1 means not set TODO define that in MetaData
  if (ival != -1) {
    if (NC::CheckErr(nc_put_att_int(ncid, varid, desc, NC_INT, 1, &ival)))
    {
      mprinterr("Error: Writing variable %s.\n", desc);
      return 1;
    }
  }
  return 0;
}

/// Add DataSet metadata to a variable
static inline int AddDataSetMetaData(MetaData const& meta, int ncid, int varid)
{
  // Filename
  if (AddDataSetStringAtt(meta.Fname().Full(), "filename", ncid, varid)) return 1;
  // Name
  if (AddDataSetStringAtt(meta.Name(), "name", ncid, varid)) return 1;
  // Aspect
  if (AddDataSetStringAtt(meta.Aspect(), "aspect", ncid, varid)) return 1;
  // Legend
  if (AddDataSetStringAtt(meta.Legend(), "legend", ncid, varid)) return 1;
  // Index
  if (AddDataSetIntAtt(meta.Idx(), "index", ncid, varid)) return 1;
  // Ensemble number
  if (AddDataSetIntAtt(meta.EnsembleNum(), "ensemblenum", ncid, varid)) return 1;
  // TODO  TimeSeries?
  // Scalar type
  if (meta.ScalarType() != MetaData::UNDEFINED) {
    if (AddDataSetStringAtt(meta.TypeString(), "scalartype", ncid, varid)) return 1;
  }
  // Scalar mode
  if (meta.ScalarMode() != MetaData::UNKNOWN_MODE) {
    if (AddDataSetStringAtt(meta.ModeString(), "scalarmode", ncid, varid)) return 1;
  }

  return 0;
}

/** Write 1D DataSets that share an index dimension. */
int DataIO_NetCDF::writeData_1D(DataSet const* ds, Dimension const& dim, SetArray const& sets) {
  mprintf("DEBUG: Sets for dimension '%s' %f %f:", dim.label(), dim.Min(), dim.Step());
  for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it)
    mprintf(" %s", it->DS()->legend());
  mprintf("\n");
  // Define the dimension. Ensure name is unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  std::string dimLabel = dim.Label() + integerToString(dimIdx_++);
  int dimId = -1;
  if (NC::CheckErr( nc_def_dim(ncid_, dimLabel.c_str(), ds->Size(), &dimId ))) {
    mprinterr("Error: Could not define dimension ID for set '%s'\n", ds->legend());
    return 1;
  }
  //int dimensionID[NC_MAX_VAR_DIMS];
  int dimensionID[1];
  dimensionID[0] = dimId;
  // Define the 'index' variable. Will be shared by this group.
  std::string idxName = dimLabel + "." + "idx";
  int idxId;
  if ( NC::CheckErr( nc_def_var(ncid_, idxName.c_str(), NC_DOUBLE, 1, dimensionID, &idxId) ) ) {
    mprinterr("Error: Could not define index variable.\n");
    return 1;
  }
  
  // Define the variable(s). Names should be unique: <DimName>.<VarName>
  for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it) {
    // Choose type
    nc_type dtype;
    switch (it->DS()->Type()) {
      case DataSet::DOUBLE  :
      case DataSet::XYMESH  : dtype = NC_DOUBLE ; break;
      case DataSet::INTEGER : dtype = NC_INT ; break;
      case DataSet::FLOAT   :
      case DataSet::PH      : dtype = NC_FLOAT ; break;
      case DataSet::UNSIGNED_INTEGER : dtype = NC_UINT ; break; // TODO netcdf4 only?
      default:
        mprinterr("Internal Error: Unhandled DataSet type for 1D NetCDF variable.\n");
        return 1;
    }
    std::string varName = dimLabel + "." + it->DS()->Meta().PrintName();
    if ( NC::CheckErr( nc_def_var(ncid_, varName.c_str(), dtype, 1, dimensionID, varIDptr_ + it->OriginalIdx() ) ) ) {
      mprinterr("Error: Could not define variable '%s'\n", varName.c_str());
      return 1;
    }
    // Add DataSet metadata as attributes
    if (AddDataSetMetaData( it->DS()->Meta(), ncid_, varIDs_[it->OriginalIdx()] )) return 1;
  } // END define variable(s)
  if (EndDefineMode( ncid_ )) return 1;
  // Write the variables
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = ds->Size();
  // Use the first set for the 'index'
  DataSet_1D const& dsidx = static_cast<DataSet_1D const&>( *(sets.front().DS()) );
  std::vector<double> idxs;
  idxs.reserve(dsidx.Size());
  for (unsigned int ii = 0; ii < dsidx.Size(); ii++)
    idxs.push_back( dsidx.Xcrd(ii) );
  if (NC::CheckErr(nc_put_vara(ncid_, idxId, start, count, &idxs[0]))) {
    mprinterr("Error: Could not write index variable from '%s'\n", dsidx.legend());
    return 1;
  }
  for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it) {
    DataSet_1D const& ds1d = static_cast<DataSet_1D const&>( *(it->DS()) );
    if (NC::CheckErr(nc_put_vara(ncid_, varIDs_[it->OriginalIdx()], start, count, ds1d.VoidPtr(0)))) {
      mprinterr("Error: Could not write variable '%s'\n", ds1d.legend());
      return 1;
    }
  }
  return 0;
}
#endif /* BINTRAJ */

// DataIO_NetCDF::WriteData()
int DataIO_NetCDF::WriteData(FileName const& fname, DataSetList const& dsl)
{
# ifdef BINTRAJ
  ncid_ = -1;
  dimIdx_ = 0;
  varIDs_.clear();
  // TODO check existing file
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid_ ) ))
    return 1;

  // Place incoming DataSets into a pool. As they are handled they will
  // be removed from the pool.
  SetPool setPool( dsl );

  varIDs_.assign( dsl.size(), -1 );
  varIDptr_ = &varIDs_[0];

  // Check our incoming data sets. Try to find common dimensions.
  for (unsigned int idx = 0; idx < setPool.Nsets(); idx++)
  {
    if (setPool.IsUsed(idx)) continue;

    DataSet const* ds = setPool.Set( idx );

    if (ds->Group() == DataSet::SCALAR_1D) {
      // ----- 1D scalar -------------------------
      SetArray sets(1, Set(ds, idx));
      setPool.MarkUsed( idx );
      // Group this set with others that share the same dimension and length.
      Dimension const& dim = ds->Dim(0);
      for (unsigned int jdx = idx + 1; jdx < setPool.Nsets(); jdx++)
      {
        DataSet const* ds2 = setPool.Set( jdx );
        if (ds->Size() == ds2->Size() && dim == ds2->Dim(0)) {
          sets.push_back( Set(ds2, jdx) );
          setPool.MarkUsed( jdx );
        }
      }
      if (writeData_1D(ds, dim, sets)) {
        mprinterr("Error: 1D NetCDF data set write failed.\n");
        return 1;
      }
    } else {
      mprinterr("Error: '%s' is an unhandled set type for NetCDF.\n", ds->legend());
      return 1;
    }
  } // END loop over set pool
  // Warn if for some reason we didnt use all the sets.
  if (!setPool.AllUsed()) {
    mprintf("Warning: Not all sets were used.\n");
  }
  mprintf("DEBUG: Variable IDs:\n");
  for (unsigned int idx = 0; idx != dsl.size(); idx++)
    mprintf("\t'%s' vid= %i\n", dsl[idx]->legend(), varIDs_[idx]);

  // Attributes
  if (EnterDefineMode(ncid_)) return 1;
  //if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"title",title.size(),title.c_str())) ) {
  //  mprinterr("Error: Writing title.\n");
  //  return 1;
  //}
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"application",5,"AMBER")) ) {
    mprinterr("Error: Writing application.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"program",7,"cpptraj")) ) {
    mprinterr("Error: Writing program.\n");
    return 1;
  }
  std::string programVersion(CPPTRAJ_INTERNAL_VERSION);
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"programVersion",
                                   programVersion.size(), programVersion.c_str())))
  {
    mprinterr("Error: Writing program version.\n");
    return 1;
  }

  // End netcdf definitions
  if (EndDefineMode(ncid_)) {
    mprinterr("NetCDF data error on ending definitions.");
    return 1;
  }
  NC::Debug(ncid_);

  nc_close(ncid_); 
  return 0;
# else
  return 1;
# endif
}
