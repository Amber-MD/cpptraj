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
// DataSets
# include "DataSet_Mesh.h"
# include "DataSet_MatrixDbl.h"
# include "DataSet_Modes.h"
# include "DataSet_3D.h"
# include "DataSet_string.h"
# include "DataSet_Vector.h"
# include "DataSet_Vector_Scalar.h"
# include "DataSet_PairwiseCache_MEM.h"
#endif

/// CONSTRUCTOR
DataIO_NetCDF::DataIO_NetCDF() :
  DataIO(true, true, true), // Valid for 1D, 2D, 3D
  ncid_(-1),
  user_specified_name_(false)
{
  SetValid( DataSet::MODES );
  SetValid( DataSet::VECTOR_SCALAR );
  SetValid( DataSet::PMATRIX_MEM );
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
    // Check that we have cpptraj conventions
    return (NC::GetConventions(infile.Filename().Full()) == NC::NC_CPPTRAJDATA);
#   else
    mprintf("Warning: '%s' is a NetCDF file but CPPTRAJ was compiled without NetCDF support.\n",
            infile.Filename().full());
#   endif
  }
  return false;
}

// -----------------------------------------------------------------------------
/** Hold info for a netcdf variable. */
class DataIO_NetCDF::NcVar {
  public:
    /// CONSTRUCTOR - blank
    NcVar() : vid_(-999), vtype_(0), hasBeenRead_(false) {}
    /// CONSTRUCTOR
    NcVar(int vidIn, nc_type vtypeIn, const char* vnameIn, int ndims, const int* dimids) :
      vid_(vidIn), vtype_(vtypeIn), vname_(vnameIn), hasBeenRead_(false), dimIds_(ndims)
    {
      for (int i = 0; i < ndims; i++) dimIds_[i] = dimids[i];
    }
    /// COPY CONSTRUCTOR
    NcVar(NcVar const& rhs) :
      vid_(rhs.vid_), vtype_(rhs.vtype_), vname_(rhs.vname_), hasBeenRead_(rhs.hasBeenRead_),
      dimIds_(rhs.dimIds_) {}
    /// ASSIGNMENT
    NcVar& operator=(NcVar const& rhs) {
      if (this == &rhs) return *this;
      vid_ = rhs.vid_;
      vtype_ = rhs.vtype_;
      vname_ = rhs.vname_;
      hasBeenRead_ = rhs.hasBeenRead_;
      dimIds_ = rhs.dimIds_;
      return *this;
    }

    int VID()                   const { return vid_; }
    nc_type Vtype()             const { return vtype_; }
    std::string const& Vname()  const { return vname_; }
    const char* vname()         const { return vname_.c_str(); }
    bool HasBeenRead()          const { return hasBeenRead_; }
    unsigned int Ndims()        const { return dimIds_.size(); }
    int DimId(unsigned int idx) const { return dimIds_[idx]; }
    bool Empty()                const { return (vid_ == -999); }
    void Print() const {
      mprintf("DEBUG:\tVariable %i - '%s' type %i dimIds:", vid_, vname_.c_str(), (int)vtype_);
      for (std::vector<int>::const_iterator it = dimIds_.begin(); it != dimIds_.end(); ++it)
        mprintf(" %i", *it);
      mprintf("\n");
    }

    void MarkRead() { hasBeenRead_ = true; }
  private:
    int vid_;                 ///< Netcdf variable id
    nc_type vtype_;           ///< Netcdf variable type
    std::string vname_;       ///< Variable name
    bool hasBeenRead_;        ///< True if the variable has been read in.
    std::vector<int> dimIds_; ///< Netcdf dimension ids for this var.
};
// -----------------------------------------------------------------------------
/** Hold info for a netcdf dimension. */
class DataIO_NetCDF::NcDim {
  public:
    /// CONSTRUCTOR - blank
    NcDim() : did_(-999), size_(0) {}
    /// CONSTRUCTOR
    NcDim(int did, std::string const& lbl, unsigned int sze) :
      did_(did), label_(lbl), size_(sze) {}
    /// COPY CONSTRUCTOR
    NcDim(NcDim const& rhs) : did_(rhs.did_), label_(rhs.label_), size_(rhs.size_) {}
    /// ASSIGNMENT
    NcDim& operator=(NcDim const& rhs) {
      if (this == &rhs) return *this;
      did_ = rhs.did_;
      label_ = rhs.label_;
      size_ = rhs.size_;
      return *this;
    }

    int DID()                  const { return did_;           }
    std::string const& Label() const { return label_;         }
    unsigned int Size()        const { return size_;          }
    bool Empty()               const { return (did_ == -999); }
    void Print() const {
      mprintf("DEBUG:\tDimension %i - '%s' (%u)\n", did_, label_.c_str(), size_);
    }
  private:
    int did_;           ///< Netcdf dimension ID
    std::string label_; ///< Dimension label
    unsigned int size_; ///< Dimension size
};

// -----------------------------------------------------------------------------

// DataIO_NetCDF::ReadHelp()
void DataIO_NetCDF::ReadHelp()
{

}

// DataIO_NetCDF::processReadArgs()
int DataIO_NetCDF::processReadArgs(ArgList& argIn)
{

  return 0;
}

/// Get integer attribute from variable
/** \return 1 if an error occurred, 0 if found, and -1 if not present. */
static inline int GetVarIntAtt(int& ival, const char* desc, int ncid, int varid)
{
  ival = -1;
  int ncerr = nc_get_att_int(ncid, varid, desc, &ival);
  if (ncerr != NC_NOERR) {
    if (ncerr != NC_ENOTATT) {
      NC::CheckErr(ncerr);
      mprinterr("Error: Could not get '%s' attribute.\n", desc);
      return 1;
    }
    return -1;
  }
  return 0;
}

/// Get double attribute from variable
static inline int GetVarDblAtt(double& dval, const char* desc, int ncid, int varid)
{
  dval = 0;
  int ncerr = nc_get_att_double(ncid, varid, desc, &dval);
  if (ncerr != NC_NOERR) {
    if (ncerr != NC_ENOTATT) {
      NC::CheckErr(ncerr);
      mprinterr("Error: Could not get '%s' attribute.\n", desc);
      return 1;
    }
    return -1;
  }
  return 0;
}

/// Get double attribute array from variable
static inline int GetVarDblArrayAtt(double* darray, const char* desc, int ncid, int varid)
{
  int ncerr = nc_get_att_double(ncid, varid, desc, darray);
  if (ncerr != NC_NOERR) {
    if (ncerr != NC_ENOTATT) {
      NC::CheckErr(ncerr);
      mprinterr("Error: Could not get '%s' attribute array.\n", desc);
      return 1;
    }
    return -1;
  }
  return 0;
}

/// Get DataSet metadata from a variable
static inline MetaData GetVarMetaData(int& errStat, int ncid, int varid)
{
  MetaData meta;
  errStat = 0;
  // Filename
  std::string att = NC::GetAttrText(ncid, varid, "filename");
  if (!att.empty())
    meta.SetFileName( att );
  // Name
  att = NC::GetAttrText(ncid, varid, "name");
  if (att.empty()) {
    mprinterr("Error: 'name' attribute missing for VID %i\n", varid);
    errStat = 1;
    return meta;
  }
  meta.SetName( att );
  // Aspect
  att = NC::GetAttrText(ncid, varid, "aspect");
  if (!att.empty())
    meta.SetAspect( att );
  // Legend
  att = NC::GetAttrText(ncid, varid, "legend");
  if (!att.empty())
    meta.SetLegend( att );
  // Index
  int ival;
  int ret = GetVarIntAtt(ival, "index", ncid, varid);
  if (ret == 1) {
    errStat = 1;
    return meta;
  } else if (ret == 0) {
    meta.SetIdx( ival );
  }
  // Ensemble number
  ret = GetVarIntAtt(ival, "ensemblenum", ncid, varid);
  if (ret == 1) {
    errStat = 1;
    return meta;
  } else if (ret == 0) {
    meta.SetEnsembleNum( ival );
  } 
  // TODO  TimeSeries?
  // Scalar mode
  MetaData::scalarMode smode = MetaData::UNKNOWN_MODE;
  att = NC::GetAttrText(ncid, varid, "scalarmode");
  if (!att.empty())
    smode = MetaData::ModeFromKeyword( att );
  meta.SetScalarMode( smode );
  // Scalar type
  MetaData::scalarType stype = MetaData::UNDEFINED;
  att = NC::GetAttrText(ncid, varid, "scalartype");
  if (!att.empty())
    stype = MetaData::TypeFromKeyword( att, smode );
  meta.SetScalarType( stype );

  return meta;
}

/// \return Variable index information from attributes
std::vector<Dimension> DataIO_NetCDF::getVarIndexInfo(int& errStat, int ncid, int varid)
const
{
//  idxVarIds.clear();
  errStat = 0;
  std::vector<Dimension> Dims;
  // Get # of index dimensions
  int nindexvar = 0;
  int ret = GetVarIntAtt(nindexvar, "nindexvar", ncid, varid);
  if (ret != 0) {
    mprinterr("Error: Missing 'nindexvar' attribute.\n");
    errStat = 1;
    return Dims;
  }
  if (debug_ > 0)
    mprintf("DEBUG: nindexvar= %i\n", nindexvar);
  if (nindexvar < 1) return Dims;

  Dims.resize( nindexvar );
  for (int i = 0; i < nindexvar; i++) {
    Dimension& dim = Dims[i];
    // Expect either label<i> and indexid<i>, or
    // label<i>, min<i>, and step<i>
    std::string suffix( integerToString(i) );
    std::string label = "label" + suffix;
    std::string att = NC::GetAttrText(ncid, varid, label.c_str());
    if (att.empty()) {
      mprintf("Warning: Index variable %i missing 'label' attribute.\n", i);
    } else {
      dim.SetLabel( att );
    }
    // Check for index<i>
    std::string str = "indexid" + suffix;
    int indexid = -1;
    int ret = GetVarIntAtt(indexid, str.c_str(), ncid, varid);
    if (ret == 1) {
      mprinterr("Error: Getting '%s' attribute.\n", str.c_str());
      errStat = 1;
      return Dims;
    } else if (ret == 0) {
//      // index<i> present
//      idxVarIds[i] = indexid;
      // TODO setting dim for backwards compatibility. This should
      // not be necessary eventually when index vars are better handled.
      dim.ChangeMin(1.0);
      dim.ChangeStep(1.0);
    } else if (ret == -1) {
      // No index<i>, check for min<i> and step<i>
      std::string min = "min" + suffix;
      std::string step = "step" + suffix;
      double dval;
      ret = GetVarDblAtt(dval, min.c_str(), ncid, varid);
      if (ret != 0) {
        mprinterr("Error: '%s' attribute is missing for varid %i.\n", min.c_str(), varid);
        errStat = 1;
        return Dims;
      }
      dim.ChangeMin( dval );
      ret = GetVarDblAtt(dval, step.c_str(), ncid, varid);
      if (ret != 0) {
        mprinterr("Error: '%s' attribute is missing for varid %i.\n", step.c_str(), varid);
        errStat = 1;
        return Dims;
      }
      dim.ChangeStep( dval );
    }
  }
  return Dims;
}

/** \return the length of the dimension of a 1D variable. */
int DataIO_NetCDF::get_1D_var_dimlen(size_t& len, int varid, const char* desc) const {
  len = 0;
  // Get the first dimension dim id
  int lengthDim[1];
  if (NC::CheckErr(nc_inq_vardimid(ncid_, varid, lengthDim))) {
    mprinterr("Error: Could not get dimension id for %s var.\n", desc);
    return 1;
  }
  // Get lengths dim size
  if (NC::CheckErr(nc_inq_dimlen(ncid_, lengthDim[0], &len))) {
    mprinterr("Error: Could not get size of dimension for %s var.\n", desc);
    return 1;
  }
  return 0;
}

/// Shortcut for getting dimension length from Dimensions_ array
unsigned int DataIO_NetCDF::dimLen(int did) const {
  return Dimensions_[did].Size();
}

/** Read CPPTRAJ XY mesh set with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_1D_xy(DataSet* ds, NcVar const& yVar, VarArray& Vars) const {
  // ----- XY Mesh ---------------
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(yVar.DimId(0));
  // Get the Y var id
  int xvarid;
  int ret = GetVarIntAtt(xvarid, "indexid0", ncid_, yVar.VID());
  if (ret != 0) {
    mprinterr("Error: No 'Yid' attribute for XY set '%s'.\n", yVar.vname());
    return 1;
  }
  DataSet_Mesh& set = static_cast<DataSet_Mesh&>( *ds );
  set.Resize(count[0]);
  DataSet_Mesh::Darray& Xvals = set.SetMeshX();
  DataSet_Mesh::Darray& Yvals = set.SetMeshY();
  // Get X
  if (NC::CheckErr(nc_get_vara(ncid_, xvarid, start, count, (void*)(&Xvals[0])))) {
    mprinterr("Error: Could not get X values for XY set.\n");
    return 1;
  }
  Vars[xvarid].MarkRead();
  // Get Y
  if (NC::CheckErr(nc_get_vara(ncid_, yVar.VID(), start, count, (void*)(&Yvals[0])))) {
    mprinterr("Error: Could not get Y values for XY set.\n");
    return 1;
  }
  Vars[yVar.VID()].MarkRead();
  return 0;
}

/** Read vector set with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_1D_vector(DataSet* ds, NcVar const& yVar, VarArray& Vars) const {
  size_t start[2];
  start[1] = 0;
  size_t count[2];
  count[0] = 1;
  count[1] = 3;

  unsigned int nvecs = dimLen(yVar.DimId(0));

  DataSet_Vector& vecSet = static_cast<DataSet_Vector&>( *ds );
  //vecSet.Allocate(DataSet::SizeArray(1, nvecs));

  // Get the xyz values
  vecSet.Resize( nvecs );
  for (unsigned int idx = 0; idx < nvecs; idx++) {
    start[0] = idx;
    if (NC::CheckErr(nc_get_vara(ncid_, yVar.VID(), start, count, vecSet[idx].Dptr()))) {
      mprinterr("Error: Could not get vector xyz values for %u '%s'\n", idx, vecSet.legend());
      return 1;
    }
  }
  Vars[yVar.VID()].MarkRead();
  // Get the origin values
  int originVarid;
  int ret = GetVarIntAtt(originVarid, "originsid", ncid_, yVar.VID());
  if (ret == 1) {
    mprinterr("Error: Problem getting 'originsid' attribute for '%s'\n", vecSet.legend());
    return 1;
  } else if (ret == 0) {
    vecSet.ResizeOrigins( nvecs );
    for (unsigned int idx = 0; idx < nvecs; idx++) {
      start[0] = idx;
      if (NC::CheckErr(nc_get_vara(ncid_, originVarid, start, count, vecSet.ModifyOxyz(idx).Dptr()))) {
        mprinterr("Error: Could not get vector origin values for %u '%s'\n", idx, vecSet.legend());
        return 1;
      }
    }
    Vars[originVarid].MarkRead();
  }
  return 0;
}

/** Read vector/scalar set with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_1D_vector_scalar(DataSet* ds, NcVar const& yVar, VarArray& Vars) const {
  size_t start[2];
  start[1] = 0;
  size_t count[2];
  count[0] = 1;
  count[1] = 3;

  unsigned int nvecs = dimLen(yVar.DimId(0));

  DataSet_Vector_Scalar& vecSet = static_cast<DataSet_Vector_Scalar&>( *ds );
  //vecSet.Allocate(DataSet::SizeArray(1, nvecs));

  // Get the xyz values
  vecSet.Resize( nvecs );
  for (unsigned int idx = 0; idx < nvecs; idx++) {
    start[0] = idx;
    if (NC::CheckErr(nc_get_vara(ncid_, yVar.VID(), start, count, vecSet.ModifyVec(idx).Dptr()))) {
      mprinterr("Error: Could not get vector xyz values for %u '%s'\n", idx, vecSet.legend());
      return 1;
    }
  }
  Vars[yVar.VID()].MarkRead();

  // Get the scalar values
  int scalarVid;
  int ret = GetVarIntAtt(scalarVid, "scalarid", ncid_, yVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'scalarid' attribute for '%s'\n", vecSet.legend());
    return 1;
  }

  start[0] = 0;
  count[0] = nvecs;
  if (NC::CheckErr(nc_get_vara(ncid_, scalarVid, start, count, vecSet.ValPtr()))) {
    mprinterr("Error: Could not get scalar values for '%s'\n", vecSet.legend());
    return 1;
  }
  Vars[scalarVid].MarkRead();

  return 0;
}

/** Read string set with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_1D_string(DataSet* ds, NcVar const& yVar, VarArray& Vars) const {
  // Get the lengths var id
  int lengthsVarid;
  int ret = GetVarIntAtt(lengthsVarid, "strlengthsid", ncid_, yVar.VID());
  if (ret != 0) {
    mprinterr("Error: No 'strlengthsid' attribute for string set '%s'.\n", yVar.vname());
    return 1;
  }
  // Get the string lengths array dimension size
  size_t nstrings;
  if (get_1D_var_dimlen(nstrings, lengthsVarid, "string lengths")) return 1;

  DataSet_string& strSet = static_cast<DataSet_string&>( *ds );
  strSet.Resize( nstrings );

  size_t cstart[1], nstart[1];
  size_t ccount[1], ncount[1];
  cstart[0] = 0;
  nstart[0] = 0;
  ncount[0] = 1;
  std::vector<char> buffer;
  for (size_t idx = 0; idx < nstrings; idx++) {
    //if (debug_ > 1) mprintf("DEBUG: read string %zu", idx);
    // Get the string length
    int len;
    if (NC::CheckErr(nc_get_vara_int(ncid_, lengthsVarid, nstart, ncount, &len))) {
      mprinterr("Error: Could not get length of string %zu\n", idx);
      return 1;
    }
    //if (debug_ > 1) mprintf(" len=%i", len);
    if ((size_t)len+1 > buffer.size())
      buffer.resize( len+1 );
    nstart[0]++;
    // Get the string
    ccount[0] = len;
    if (NC::CheckErr(nc_get_vara_text(ncid_, yVar.VID(), cstart, ccount, &buffer[0]))) {
      mprinterr("Error: Could not get string %zu\n", idx);
      return 1;
    }
    buffer[len] = '\0';
    //if (debug_ > 1) mprintf(" str='%s'\n", &buffer[0]);
    cstart[0] += (size_t)len;
    // Assign to DataSet
    strSet.SetElement( idx, &buffer[0] );
  }
  // Mark vars as read
  Vars[yVar.VID()].MarkRead();
  Vars[lengthsVarid].MarkRead();

  return 0;
}

/** Read 1D array with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_1D(DataSet* ds, NcVar const& yVar, VarArray& Vars) const {
  // ----- 1D Scalar -------------
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(yVar.DimId(0));

  DataSet_1D& set = static_cast<DataSet_1D&>( *ds );
  set.Resize( count[0] );
  if (NC::CheckErr(nc_get_vara(ncid_, yVar.VID(), start, count, (void*)(set.Yptr())))) {
    mprinterr("Error: Could not get values for set.\n");
    return 1;
  }
  Vars[yVar.VID()].MarkRead();
  return 0;
}

/** Read cluster pairwise matrix with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_cluster_pwmatrix(DataSet* ds, NcVar const& matrixVar, VarArray& Vars) const {
  // Get n_original_frames attribute
  int n_original_frames;
  int ret = GetVarIntAtt(n_original_frames, "n_original_frames", ncid_, matrixVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'n_original_frames' attribute.\n");
    return 1;
  }
  // Get metric description
  std::string metric_description = NC::GetAttrText(ncid_, matrixVar.VID(), "MetricDescription");
  // Get the sieve attribute
  int sieve;
  ret = GetVarIntAtt(sieve, "sieve", ncid_, matrixVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'sieve' attribute.\n");
    return 1;
  }
  // Get the frames being clustered
  Cpptraj::Cluster::Cframes frames_to_cluster; 
  int actual_framesid = -1;
  if (sieve != 1) {
    ret = GetVarIntAtt(actual_framesid, "indexid0", ncid_, matrixVar.VID());
    if (ret != 0) {
      mprinterr("Error: Could not get 'indexid0' attribute (for clustered frames variable ID) but sieve != 1.\n");
      return 1;
    }
    // TODO check that indexid0 == indexid1?
    // Get the actual_frames dimension length (i.e. # rows) TODO check against nrows attribute?
    size_t n_actual_frames;
    if (get_1D_var_dimlen(n_actual_frames, actual_framesid, "actual # frames")) return 1;
    // Read actual frames
    frames_to_cluster.assign(n_actual_frames, -1);
    size_t start[1], count[1];
    start[0] = 0;
    count[0] = n_actual_frames;
    if (NC::CheckErr(nc_get_vara(ncid_, actual_framesid, start, count, frames_to_cluster.Ptr()))) {
      mprinterr("Error: Could not read actual frames array.\n");
      return 1;
    }
    Vars[actual_framesid].MarkRead();
  } else {
    // Sieve 1 (i.e. no sieveing)
    frames_to_cluster.reserve(n_original_frames);
    for (int frm = 0; frm < n_original_frames; frm++)
      frames_to_cluster.push_back( frm );
  }
  if (debug_ > 1) {
    mprintf("DEBUG: Frames in pairwise matrix:");
    for (Cpptraj::Cluster::Cframes::const_iterator it = frames_to_cluster.begin();
                                                   it != frames_to_cluster.end(); ++it)
      mprintf(" %i", *it);
    mprintf("\n");
  }
  // Allocate matrix
  DataSet_PairwiseCache_MEM& pwmatrix = static_cast<DataSet_PairwiseCache_MEM&>( *ds );
  if (pwmatrix.SetupCache( n_original_frames, frames_to_cluster, sieve, metric_description ))
  {
    mprinterr("Error: Failed to set up cluster pairwise matrix.\n");
    return 1;
  }
  // Read matrix
  size_t start[1], count[1];
  start[0] = 0;
  count[0] = dimLen(matrixVar.DimId(0));
  if (NC::CheckErr(nc_get_vara(ncid_, matrixVar.VID(), start, count, pwmatrix.Ptr()))) {
    mprinterr("Error: Could not read pairwise matrix variable.\n");
    return 1;
  }

  return 0;
}

/** Read 2D matrix with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_2D(DataSet* ds, NcVar const& matVar, VarArray& Vars) const {
  // ----- 2D Matrix -------------
  // Get nrows/ncols
  int ncols, nrows;
  int ret = GetVarIntAtt(ncols, "ncols", ncid_, matVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'ncols'.\n");
    return 1;
  }
  ret = GetVarIntAtt(nrows, "nrows", ncid_, matVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'nrows'.\n");
    return 1;
  }
  // Get matrix kind
  std::string mkind = NC::GetAttrText(ncid_, matVar.VID(), "matrixkind");
  // Allocate
  DataSet_2D& mat = static_cast<DataSet_2D&>( *ds );
  int allocErr = 0;
  if (mkind == "full")
    allocErr = mat.Allocate2D(ncols, nrows);
  else if (mkind == "half")
    allocErr = mat.AllocateHalf(ncols);
  else if (mkind == "tri")
    allocErr = mat.AllocateTriangle(ncols);
  else {
    mprinterr("Error: Urecognized matrix kind: %s\n", mkind.c_str());
    return 1;
  }
  if (allocErr != 0) {
    mprinterr("Error: Could not allocate matrix.\n");
    return 1;
  }
  // Read values
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(matVar.DimId(0));
  if (NC::CheckErr(nc_get_vara(ncid_, matVar.VID(), start, count, (void*)(mat.MatrixPtr())))) {
    mprinterr("Error: Could not get values for matrix.\n");
    return 1;
  }
  Vars[matVar.VID()].MarkRead();
  // Check for nsnapshots
  int nsnapshots = 0;
  ret = GetVarIntAtt(nsnapshots, "nsnapshots", ncid_, matVar.VID());
  if (ret == 1)
    return 1;
  else if (ret == 0) {
    DataSet_MatrixDbl& dmat = static_cast<DataSet_MatrixDbl&>( mat );
    dmat.SetNsnapshots( nsnapshots );
  }
  // Check for vectid
  int vectVarId = -1;
  ret = GetVarIntAtt(vectVarId, "vectid", ncid_, matVar.VID());
  if (ret == 1) return 1;
  if (ret == 0) {
    mprintf("\tMatrix has diagonal vector data.\n");
    if (mat.Type() != DataSet::MATRIX_DBL) {
      mprinterr("Error: Variable has vect data but set is not double matrix.\n");
      return 1;
    }
    unsigned int vectLength = dimLen(Vars[vectVarId].DimId(0));
    DataSet_MatrixDbl& dmat = static_cast<DataSet_MatrixDbl&>( mat );
    dmat.AllocateVector( vectLength );
    count[0] = vectLength;
    if (NC::CheckErr(nc_get_vara(ncid_, vectVarId, start, count, (void*)(&dmat.V1()[0])))) {
      mprinterr("Error: Could not get vect for matrix.\n");
      return 1;
    }
    Vars[vectVarId].MarkRead();
  }
  // Check for massid
  int massVarId = -1;
  ret = GetVarIntAtt(massVarId, "massid", ncid_, matVar.VID());
  if (ret == 1) return 1;
  if (ret == 0) {
    mprintf("\tMatrix has mass data.\n");
    if (mat.Type() != DataSet::MATRIX_DBL) {
      mprinterr("Error: Variable has mass data but set is not double matrix.\n");
      return 1;
    }
    unsigned int massLength = dimLen(Vars[massVarId].DimId(0));
    DataSet_MatrixDbl& dmat = static_cast<DataSet_MatrixDbl&>( mat );
    dmat.AllocateMass( massLength );
    count[0] = massLength;
    if (NC::CheckErr(nc_get_vara(ncid_, massVarId, start, count, (void*)(&dmat.M1()[0])))) {
      mprinterr("Error: Could not get mass for matrix.\n");
      return 1;
    }
    Vars[massVarId].MarkRead();
  }
  return 0;
}

/** Read 3D grid with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_3D(DataSet* ds, NcVar const& gridVar, VarArray& Vars) const {
  // ----- 3D Grid ---------------------
  // Get nx/ny/nz 
  int nx, ny, nz;
  int ret = GetVarIntAtt(nx, "nx", ncid_, gridVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'nx'.\n");
    return 1;
  }
  ret = GetVarIntAtt(ny, "ny", ncid_, gridVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'ny'.\n");
    return 1;
  }
  ret = GetVarIntAtt(nz, "nz", ncid_, gridVar.VID());
  if (ret != 0) {
    mprinterr("Error: Could not get 'nz'.\n");
    return 1;
  }
  // Get origin
  double oxyz[3];
  ret = GetVarDblArrayAtt(oxyz, "origin", ncid_, gridVar.VID());
  if (ret == 1) {
    mprinterr("Error: Could not get attribute 'origin'.\n");
    return 1;
  } else if (ret == -1) {
    mprinterr("Error: Missing attribute 'origin'.\n");
    return 1;
  }
  // Get ucell
  double ucell[9];
  ret = GetVarDblArrayAtt(ucell, "ucell", ncid_, gridVar.VID());
  if (ret == 1) {
    mprinterr("Error: Could not get attribute 'ucell'.\n");
    return 1;
  } else if (ret == -1) {
    mprinterr("Error: Missing attribute 'ucell'.\n");
    return 1;
  }
  Box gridBox;
  if (gridBox.SetupFromUcell( ucell )) {
    mprinterr("Error: Could not set up grid unit cell.\n");
    return 1;
  }

  // Allocate
  DataSet_3D& grid = static_cast<DataSet_3D&>( *ds );
  int allocErr = grid.Allocate_N_O_Box( nx, ny, nz, Vec3(oxyz), gridBox );
  if (allocErr != 0) {
    mprinterr("Error: Could not allocate grid.\n");
    return 1;
  }
  // Read values
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(gridVar.DimId(0));
  if (NC::CheckErr(nc_get_vara(ncid_, gridVar.VID(), start, count, (void*)(grid.GridPtr())))) {
    mprinterr("Error: Could not get values for grid.\n");
    return 1;
  }
  Vars[gridVar.VID()].MarkRead();
  // DEBUG
  if (debug_ > 0)
    grid.GridInfo();

  return 0;
}

/** Read modes data with CPPTRAJ conventions. */
int DataIO_NetCDF::readData_modes(DataSet* ds, NcVar const& modesVar, VarArray& Vars) const {
  // ----- Modes -----------------
  unsigned int n_eigenvalues = dimLen(modesVar.DimId(0));
  // Get the eigenvectors variable ID
  int vectorsVarId = -1;
  int ret = GetVarIntAtt(vectorsVarId, "eigenvectorsid", ncid_, modesVar.VID());
  if (ret != 0) {
    mprinterr("Error: Modes data missing 'eigenvectorsid'.\n");
    return 1;
  }
  unsigned int evectorLength = dimLen(Vars[vectorsVarId].DimId(0));
  // Get the avg. coords variable ID
  int coordsVarId = -1;
  ret = GetVarIntAtt(coordsVarId, "avgcoordsid", ncid_, modesVar.VID());
  if (ret != 0) {
    mprinterr("Error: Modes data missing 'avgcoordsid'.\n");
    return 1;
  }
  unsigned int avgCoordsLength = dimLen(Vars[coordsVarId].DimId(0));
  // Get the mass variable ID if present
  int massVarId = -1;
  unsigned int massLength = 0;
  ret = GetVarIntAtt(massVarId, "massid", ncid_, modesVar.VID());
  if (ret == 1) {
    mprinterr("Error: Could not get 'massid'.\n");
    return 1;
  } else if (ret == 0) {
    massLength = dimLen(Vars[massVarId].DimId(0));
  }
  mprintf("\tModes: # values= %u, evector size= %u, avg coords size= %u, mass size = %u\n",
          n_eigenvalues, evectorLength, avgCoordsLength, massLength);
  // Allocate the modes set
  DataSet_Modes& modes = static_cast<DataSet_Modes&>( *ds );
  if (modes.AllocateModes(n_eigenvalues, evectorLength, avgCoordsLength, massLength)) {
    mprinterr("Error: Could not allocate memory for modes set.\n");
    return 1;
  }
  // Read modes data
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(modesVar.DimId(0));
  if (NC::CheckErr(nc_get_vara(ncid_, modesVar.VID(), start, count, modes.EvalPtr()))) {
    mprinterr("Error: Could not read eigenvalues.\n");
    return 1;
  }
  count[0] = evectorLength;
  if (NC::CheckErr(nc_get_vara(ncid_, vectorsVarId, start, count, modes.EvectPtr()))) {
    mprinterr("Error: Could not read eigenvectors.\n");
    return 1;
  }
  count[0] = avgCoordsLength;
  if (NC::CheckErr(nc_get_vara(ncid_, coordsVarId, start, count, modes.AvgFramePtr()))) {
    mprinterr("Error: Could not read avg. coords.\n");
    return 1;
  }
  if (massVarId != -1 && massLength > 0) {
    count[0] = massLength;
    if (NC::CheckErr(nc_get_vara(ncid_, massVarId, start, count, modes.MassPtr()))) {
      mprinterr("Error: Could not read masses.\n");
      return 1;
    }
  }
  // Mark variables as read
  Vars[modesVar.VID()].MarkRead();
  Vars[vectorsVarId].MarkRead();
  Vars[coordsVarId].MarkRead();
  if (massVarId != -1)
    Vars[massVarId].MarkRead();
  return 0;
}

/** Read variable with CPPTRAJ conventions. */
int DataIO_NetCDF::read_cpptraj_vars(DataSetList& dsl, std::string const& dsname, VarArray& Vars)
const
{
  for (VarArray::const_iterator var = Vars.begin(); var != Vars.end(); ++var)
  {
    if (var->HasBeenRead()) continue;

    // Get the description
    std::string desc = NC::GetAttrText(ncid_, var->VID(), "description");
    // Get the type from the description
    DataSet::DataType dtype = DataSet::TypeFromDescription( desc );
    if (debug_ > 0)
      mprintf("\t%s Description: %s (%i)\n", var->vname(), desc.c_str(), (int)dtype);
    // Get metaData
    int errStat = 0;
    MetaData meta = GetVarMetaData( errStat, ncid_, var->VID() );
    if (errStat != 0) {
      mprinterr("Error: Could not set up meta data for variable '%s'\n", var->vname());
      return 1;
    }
    // Get DataSet dimensions
    errStat = 0;
    std::vector<Dimension> Dims = getVarIndexInfo(errStat, ncid_, var->VID());
    if (errStat != 0) return 1;
    if (debug_ > 0)
      for (std::vector<Dimension>::const_iterator dim = Dims.begin(); dim != Dims.end(); ++dim)
        mprintf("DEBUG:\t Var %s dim %s min %f step %f\n", var->vname(), dim->label(), dim->Min(), dim->Step());
    // Add the set
    DataSet* ds = dsl.AddSet( dtype, meta );
    if (ds == 0) {
      mprinterr("Error: Could not allocate set '%s'\n", meta.PrintName().c_str());
      return 1;
    }
    // Set DataSet dimensions
    unsigned int idx = 0;
    for (std::vector<Dimension>::const_iterator dim = Dims.begin(); dim != Dims.end(); ++dim)
      ds->SetDim(idx++, *dim);
    // Check netcdf variable dimensions
//    if (var->Ndims() == 1) {
      // One flat dimension
      if (debug_ > 0)
        mprintf("DEBUG: %s dim length %u\n", var->vname(), dimLen(var->DimId(0)) );
      if (dtype == DataSet::XYMESH) {
        if (readData_1D_xy(ds, *var, Vars))
          return 1;
      } else if (dtype == DataSet::VECTOR) {
        if (readData_1D_vector(ds, *var, Vars))
          return 1;
      } else if (dtype == DataSet::VECTOR_SCALAR) {
        if (readData_1D_vector_scalar(ds, *var, Vars))
          return 1;
      } else if (dtype == DataSet::STRING) {
        if (readData_1D_string(ds, *var, Vars))
          return 1;
      } else if (ds->Group() == DataSet::SCALAR_1D) {
        if (readData_1D(ds, *var, Vars))
          return 1;
      } else if (dtype == DataSet::PMATRIX_MEM) {
        if (readData_cluster_pwmatrix(ds, *var, Vars))
         return 1;
      } else if (ds->Group() == DataSet::MATRIX_2D) {
        if (readData_2D(ds, *var, Vars))
          return 1;
      } else if (ds->Group() == DataSet::GRID_3D) {
        if (readData_3D(ds, *var, Vars))
          return 1;
      } else if ( dtype == DataSet::MODES ) {
        if (readData_modes(ds, *var, Vars))
          return 1;
      } else {
        mprinterr("Error: Cannot read type '%s' yet.\n", desc.c_str());
        return 1;
      }
//    } else {
//      mprinterr("Error: Cannot read type '%s' yet.\n", desc.c_str());
//      return 1;
//    }
  }
  return 0;
}

// DataIO_NetCDF::ReadData()
int DataIO_NetCDF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
# ifdef BINTRAJ
  // Check if the user specified a data set name. The default is to use the
  // file base name, optionally plus '_<index>'
  if ( dsname.find(fname.Base()) == std::string::npos ) {
    mprintf("\tUser has specified a data set name.\n");
    user_specified_name_ = true;
  } else {
    mprintf("\tUser has not specified a data set name.\n");
    user_specified_name_ = false;
  }

  if (NC::CheckErr( nc_open( fname.full(), NC_NOWRITE, &ncid_ ) != NC_NOERR )) {
    mprinterr("Error: Could not open NetCDF data file '%s'\n", fname.full());
    return 1;
  }

  // Check if we have CPPTRAJ conventions
  bool hasCpptrajConventions = (NC::GetConventions(ncid_) == NC::NC_CPPTRAJDATA);
  if (hasCpptrajConventions)
    mprintf("\tNetCDF data file has CPPTRAJ conventions.\n");
  else {
    mprinterr("Error: NetCDF data file has unknown conventions.\n");
    return 1;
  }

  // Get the number of dimensions, variables, attributes, and ID of the
  // unlimited dimension (if any).
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  if (NC::CheckErr( nc_inq(ncid_, &ndimsp, &nvarsp, &ngattsp, &unlimdimidp) )) {
    mprinterr("Error: Could not get NetCDF data file information.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: '%s' : ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            fname.full(), ndimsp, nvarsp, ngattsp, unlimdimidp);
  char varName[NC_MAX_NAME+1];
  // Get the length of all dimensions
  Dimensions_.clear();
  Dimensions_.reserve(ndimsp);
  for (int idim = 0; idim < ndimsp; idim++) {
    size_t diml;
    if (NC::CheckErr(nc_inq_dim(ncid_, idim, varName, &diml))) {
      mprinterr("Error: Could not get length of NetCDF data dimension %i\n", idim);
      return 1;
    }
    Dimensions_.push_back( NcDim( idim, varName, diml ) );
    if (debug_ > 0)
      Dimensions_.back().Print();
  }

  VarArray AllVars;
  AllVars.reserve( nvarsp ); 

  // Loop over all variables in the NetCDF file.
  nc_type varType = 0;            // Variable type
  int nVarDims = -1;              // # variable dimensions
  int varDimIds[NC_MAX_VAR_DIMS]; // variable dim ids
  int nVarAttributes = -1;        // number of variable attributes
  for (int ivar = 0; ivar < nvarsp; ivar++) {
    if (NC::CheckErr(nc_inq_var(ncid_, ivar, varName, &varType, &nVarDims, varDimIds, &nVarAttributes))) {
      mprinterr("Error: Could not get NetCDF data variable name %i\n", ivar);
      return 1;
    }
    if (debug_ > 0)
      mprintf("DEBUG:\tVariable %i - '%s', %i dims, %i attributes\n", ivar, varName, nVarDims, nVarAttributes);
    AllVars.push_back( NcVar(ivar, varType, varName, nVarDims, varDimIds) );
  }

  for (VarArray::const_iterator it = AllVars.begin(); it != AllVars.end(); ++it)
    mprintf("\t  %i (%s)\n", it->VID(), it->vname());
  
  if (read_cpptraj_vars(dsl, dsname, AllVars)) return 1;

  nc_close( ncid_ );
  return 0;
# else  
  return 1;
# endif
}

// =============================================================================
// DataIO_NetCDF::WriteHelp()
void DataIO_NetCDF::WriteHelp()
{

}

// DataIO_NetCDF::processWriteArgs()
int DataIO_NetCDF::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// -----------------------------------------------------------------------------
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
    /// Print unused sets to stdout
    void PrintUnused() const {
      for (unsigned int idx = 0; idx < sets_.size(); idx++)
        if (!isUsed_[idx]) mprintf("\tUnused: %s\n", sets_[idx]->legend());
    }
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
// -----------------------------------------------------------------------------

/** Hold a pointer to DataSet and its original index. Used to refer back to
  * original DataSetList from sets in a SetPool.
  */
class DataIO_NetCDF::Set {
  public:
    //Set(DataSet const* ds, int oidx) : ds_(ds), oidx_(oidx) {}
    Set(DataSet const* ds) : ds_(ds) {}

    DataSet const* DS() const { return ds_; }
    //int OriginalIdx() const { return oidx_; }
  private:
    DataSet const* ds_;
    //int oidx_;
};
// -----------------------------------------------------------------------------

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
      mprinterr("Error: Writing attribute %s.\n", desc);
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
      mprinterr("Error: Writing attribute %s.\n", desc);
      return 1;
    }
  }
  return 0;
}

/// Add a double as double attribute
static inline int AddDataSetDblAtt(double dval, const char* desc, int ncid, int varid)
{
  if (NC::CheckErr(nc_put_att_double(ncid, varid, desc, NC_DOUBLE, 1, &dval)))
  {
    mprinterr("Error: Writing attribute %s.\n", desc);
    return 1;
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

/// Add DataSet index information to a target variable
/** \param ds Input DataSet
  * \param indexVarIds Optional index variable IDs, one for each dimension.
  * \param ncid File ncid
  * \param varid Target variable id
  */
static inline int AddDataSetIndexInfo(DataSet const* ds, std::vector<int> const& indexVarIds,
                                      int ncid, int varid)
{
  // Add number of index variables 
  if (AddDataSetIntAtt( ds->Ndim(), "nindexvar", ncid, varid )) return 1;
  // Loop over index variables
  for (int idx = 0 ; idx < (int)ds->Ndim(); idx++)
  {
    Dimension const& dim = ds->Dim(idx);
    std::string suffix( integerToString(idx) );
    if (!indexVarIds.empty() && indexVarIds[idx] > -1) {
      // Index variable ID present for this dimension
      std::string indexid( "indexid" + suffix );
      if (AddDataSetIntAtt(indexVarIds[idx], indexid.c_str(), ncid, varid)) return 1;
    } else {
      // No index variable ID for this dimension
      std::string min(  "min"   + suffix);
      std::string step( "step"  + suffix);
      if (AddDataSetDblAtt(dim.Min(),  min.c_str(),  ncid, varid)) return 1;
      if (AddDataSetDblAtt(dim.Step(), step.c_str(), ncid, varid)) return 1;
    }
    std::string label("label" + suffix);
    if (AddDataSetStringAtt(dim.Label(), label.c_str(), ncid, varid)) return 1;
  }
  return 0;
}

/// Add DataSet MetaData, index information (with optional index var IDs), and description to target variable
/** \param ds Set with MetaData to add.
  * \param indexVarIds Optional index variable IDs, one for each dimension.
  * \param ncid NetCDF id to add to.
  * \param varid NetCDF variable id to add to.
  */
static inline int AddDataSetInfo(DataSet const* ds, std::vector<int> const& indexVarIds,
                                 int ncid, int varid)
{
  // Add DataSet metadata as attributes
  if (AddDataSetMetaData(ds->Meta(), ncid, varid)) return 1;
  // Add index info
  if (AddDataSetIndexInfo(ds, indexVarIds, ncid, varid)) return 1;
  // Store the description
  if (AddDataSetStringAtt(ds->description(), "description", ncid, varid)) return 1;
  return 0;
}

/// Add DataSet MetaData, index information (no index var IDs), and description to target variable
/** \param ds Set with MetaData to add.
  * \param ncid NetCDF id to add to.
  * \param varid NetCDF variable id to add to.
  */
static inline int AddDataSetInfo(DataSet const* ds, int ncid, int varid) {
  return AddDataSetInfo(ds, std::vector<int>(), ncid, varid);
}

/** Define dimension. Ensure name is unique by appending an index.
  * Dimension label will be '<label>.<index>'.
  * \return index of defined dimension in Dimensions_, -1 on error.
  */
int DataIO_NetCDF::defineDim(std::string const& label, unsigned int dimSize, std::string const& setname)
{
  if (label.empty()) {
    mprinterr("Internal Error: defineDim(): label is empty.\n");
    return -1;
  }
  int dimIdx = Dimensions_.size();
  std::string dimLabel( label + "." + integerToString(dimIdx) );
  int did = -1;
  if (NC::CheckErr( nc_def_dim(ncid_, dimLabel.c_str(), dimSize, &did ))) {
    mprinterr("Error: Could not define dimension '%s' ID for set '%s' (%s)\n", label.c_str(), setname.c_str());
    return -1;
  }
  Dimensions_.push_back( NcDim(did, dimLabel, dimSize) );
  if (debug_ > 0)
    Dimensions_.back().Print();
  return dimIdx;
}

/** Create ND variable. Optionally associate it with a parent variable.
  * \return Created variable.
  */
DataIO_NetCDF::NcVar DataIO_NetCDF::defineVar(std::vector<int> const& DimIds,
                                              int nctype,
                                              std::string const& PrintName,
                                              std::string const& VarSuffix,
                                              int parentVarId)
const
{
  // ASSUME WE ARE IN DEFINE MODE
  std::string varName = PrintName + "." + VarSuffix;
  int varid;
  if ( NC::CheckErr(nc_def_var(ncid_, varName.c_str(), nctype, DimIds.size(), &DimIds[0], &varid)) ) {
    mprinterr("Error: Could not define variable '%s'\n", varName.c_str());
    return NcVar();
  }
  if (parentVarId > -1) {
    std::string childVarDesc = VarSuffix + "id";
    if (AddDataSetIntAtt( varid, childVarDesc.c_str(), ncid_, parentVarId )) return NcVar();
  }
  
  NcVar ret(varid, nctype, varName.c_str(), DimIds.size(), &DimIds[0]);
  if (debug_ > 0)
    ret.Print(); // DEBUG
  return ret;
}

/** Create 1D variable. Optionally associate it with a parent variable.
  * \return Created variable.
  */
DataIO_NetCDF::NcVar DataIO_NetCDF::defineVar(int dimid, int nctype,
                                              std::string const& PrintName,
                                              std::string const& VarSuffix,
                                              int parentVarId)
const
{
  std::vector<int> DimIds(1, dimid);
  return defineVar(DimIds, nctype, PrintName, VarSuffix, parentVarId);
}

/** Create 1D variable. */
DataIO_NetCDF::NcVar DataIO_NetCDF::defineVar(int dimid, int nctype, std::string const& PrintName, std::string const& VarSuffix)
const
{
  return defineVar(dimid, nctype, PrintName, VarSuffix, -1);
}

/** Write 1D VECTOR_SCALAR set. TODO should be covered in writeData_1D. */
int DataIO_NetCDF::writeData_vector_scalar(DataSet const* ds) {
  if (EnterDefineMode(ncid_)) return 1;
  // Define the dimension. Ensure name is unique by appending an index.
  int dimIdx = defineDim( "length", ds->Size(), ds->Meta().Legend() );
  if (dimIdx < 0) return 1;
  NcDim lengthDim = Dimensions_[dimIdx];

  DataSet_Vector_Scalar const& vecSet = static_cast<DataSet_Vector_Scalar const&>( *ds );
  // Define xyz dimension
  dimIdx = defineDim("xyz", 3, vecSet.Meta().Legend() + " xyz" );
  if (dimIdx < 0) return 1;
  NcDim xyzDim = Dimensions_[dimIdx];
  // Define vector variable
  std::vector<int> DimIds(2);
  DimIds[0] = lengthDim.DID();
  DimIds[1] = xyzDim.DID();
  NcVar vecVar = defineVar( DimIds, NC_DOUBLE, vecSet.Meta().PrintName(), "vecs", -1);
  if (vecVar.Empty()) return 1;
  if (AddDataSetInfo( ds, ncid_, vecVar.VID() )) return 1;
  // Define scalar variable
  NcVar scalarVar = defineVar( lengthDim.DID(), NC_DOUBLE, vecSet.Meta().PrintName(), "scalar", vecVar.VID() );
  if (scalarVar.Empty()) return 1;
  if (EndDefineMode( ncid_ )) return 1;

  size_t start[2], count[2];
  start[1] = 0;
  count[0] = 1;
  count[1] = 3;
  for (unsigned int idx = 0; idx < vecSet.Size(); idx++) {
    start[0] = idx;
    if (NC::CheckErr(nc_put_vara(ncid_, vecVar.VID(), start, count, vecSet.Vec(idx).Dptr()))) {
      mprinterr("Error: Could not write vector %u xyz from '%s'\n", idx, vecSet.legend());
      return 1;
    }
  }

  start[0] = 0;
  count[0] = lengthDim.Size();
  if (NC::CheckErr(nc_put_vara(ncid_, scalarVar.VID(), start, count, vecSet.ValPtr()))) {
    mprinterr("Error: Could not write scalar values from '%s'\n", vecSet.legend());
    return 1;
  }

  return 0;
}

/** Write 1D DataSets that share an index dimension. */
int DataIO_NetCDF::writeData_1D(DataSet const* ds, Dimension const& dim, SetArray const& sets) {
  if (ds->Type() == DataSet::PH)
    mprintf("Warning: Currently only State information saved for pH sets.\n");
  if (debug_ > 0) {
    mprintf("DEBUG: Sets for dimension '%s' %f %f:", dim.label(), dim.Min(), dim.Step());
    for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it)
      mprintf(" %s", it->DS()->legend());
    mprintf("\n");
  }
  // Define the dimension. Ensure name is unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  int dimIdx = defineDim( "length", ds->Size(), ds->Meta().Legend() );
  if (dimIdx < 0) return 1;
  NcDim lengthDim = Dimensions_[dimIdx];
  
  // Define the variable(s). Names should be unique. May be more than
  // one variable per DataSet.
  std::vector<VarArray> variables;
  for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it) {
    // Choose type
    nc_type dtype;
    switch (it->DS()->Type()) {
      case DataSet::VECTOR_SCALAR :
      case DataSet::VECTOR  :
      case DataSet::DOUBLE  :
      case DataSet::XYMESH  :          dtype = NC_DOUBLE ; break;
      case DataSet::PH      :
      case DataSet::INTEGER :          dtype = NC_INT ; break;
      case DataSet::FLOAT   :          dtype = NC_FLOAT ; break;
      case DataSet::UNSIGNED_INTEGER : dtype = NC_UINT ; break; // TODO netcdf4 only?
      case DataSet::STRING  :          dtype = NC_CHAR ; break;
      default:
        mprinterr("Internal Error: '%s': Unhandled DataSet type for 1D NetCDF variable.\n", it->DS()->legend());
        return 1;
    }
    // Add variable(s)
    variables.push_back( VarArray() );
    VarArray& set_vars = variables.back();
    if (it->DS()->Type() == DataSet::STRING) {
      // ----- String set --------------
      DataSet_string const& strSet = static_cast<DataSet_string const&>( *(it->DS()) );
      // Sum up all characters for the netcdf character array dimension
      unsigned int nchars = 0;
      for (unsigned int idx = 0; idx < strSet.Size(); idx++)
        nchars += strSet[idx].size();
      dimIdx = defineDim("nchars", nchars, strSet.Meta().Legend() + " nchars" );
      if (dimIdx < 0) return 1;
      NcDim ncharsDim = Dimensions_[dimIdx];
      // Define the character array variable
      NcVar charsVar = defineVar(ncharsDim.DID(), dtype, strSet.Meta().PrintName(), "chars");
      if (charsVar.Empty()) return 1;
      if (AddDataSetInfo( it->DS(), ncid_, charsVar.VID() )) return 1;
      set_vars.push_back( charsVar );
      // Define the string lengths variable
      NcVar lengthsVar = defineVar(lengthDim.DID(), NC_INT, strSet.Meta().PrintName(), "strlengths", charsVar.VID());
      if (lengthsVar.Empty()) return 1;
      set_vars.push_back( lengthsVar );
    } else if (it->DS()->Type() == DataSet::VECTOR) {
      // ----- Vector set --------------
      DataSet_Vector const& vecSet = static_cast<DataSet_Vector const&>( *(it->DS()) );
      // Define xyz dimension
      dimIdx = defineDim("xyz", 3, vecSet.Meta().Legend() + " xyz" );
      if (dimIdx < 0) return 1;
      NcDim xyzDim = Dimensions_[dimIdx];
      // Define vector variable
      std::vector<int> DimIds(2);
      DimIds[0] = lengthDim.DID();
      DimIds[1] = xyzDim.DID();
      NcVar vecVar = defineVar( DimIds, dtype, vecSet.Meta().PrintName(), "vecs", -1);
      if (vecVar.Empty()) return 1;
      if (AddDataSetInfo( it->DS(), ncid_, vecVar.VID() )) return 1;
      set_vars.push_back( vecVar );
      // Define the origin variable if needed
      if (vecSet.HasOrigins()) {
        NcVar originVar = defineVar( DimIds, dtype, vecSet.Meta().PrintName(), "origins", vecVar.VID());
        set_vars.push_back( originVar );
      }
    } else if (it->DS()->Type() == DataSet::XYMESH) {
      // ----- XY mesh set -------------
      // Define the Y variable
      NcVar yVar = defineVar(Dimensions_[dimIdx].DID(), NC_DOUBLE, ds->Meta().PrintName(), "Y");
      if (yVar.Empty()) {
        mprinterr("Error: Could not define Y variable for set '%s'\n", ds->legend());
        return 1;
      }
      set_vars.push_back( yVar );
      // Define the X variable.
      NcVar xVar = defineVar(Dimensions_[dimIdx].DID(), NC_DOUBLE, ds->Meta().PrintName(), "X");
      if (xVar.Empty()) {
        mprinterr("Error: Could not define Y variable for '%s'\n", ds->legend());
        return 1;
      }
      set_vars.push_back( xVar );
      // Add DataSet info to variable
      if (AddDataSetInfo(ds, std::vector<int>(1, xVar.VID()), ncid_, yVar.VID())) return 1;
    } else {
      // ----- All other 1D sets -------
      set_vars.push_back( defineVar(lengthDim.DID(), dtype, it->DS()->Meta().PrintName(), "Y") );
      if (set_vars.back().Empty()) {
        mprinterr("Error: Could not define variable for set '%s'\n", it->DS()->legend());
        return 1;
      }
      NcVar const& currentVar = set_vars.back();
      // Add DataSet info to variable
      if (AddDataSetInfo( it->DS(), ncid_, currentVar.VID() )) return 1;
    }
  } // END define variable(s)
  if (EndDefineMode( ncid_ )) return 1;
  // Write the variables
  SetArray::const_iterator dset = sets.begin();
  for (std::vector<VarArray>::const_iterator it = variables.begin(); it != variables.end(); ++it, ++dset)
  {
    if (dset->DS()->Type() == DataSet::STRING) {
      // ----- String set --------------
      VarArray const& set_vars = *it;
      DataSet_string const& strSet = static_cast<DataSet_string const&>( *(dset->DS()) );
      size_t cstart[1], nstart[1];
      size_t ccount[1], ncount[1];
      // Write the strings. To avoid allocating a lot of extra memory write
      // one string at a time.
      cstart[0] = 0;
      ncount[0] = 1;
      for (unsigned int idx = 0; idx < strSet.Size(); idx++) {
        int len = strSet[idx].size();
        ccount[0] = strSet[idx].size();
        if (NC::CheckErr(nc_put_vara(ncid_, set_vars[0].VID(), cstart, ccount, strSet[idx].c_str()))) {
          mprinterr("Error: Could not write characters for string %u from '%s'\n", idx+1, strSet.legend());
          return 1;
        }
        cstart[0] += strSet[idx].size();

        nstart[0] = idx;
        if (NC::CheckErr(nc_put_vara_int(ncid_, set_vars[1].VID(), nstart, ncount, &len))) {
          mprinterr("Error: Could not write length for string %u from '%s'\n", idx+1, strSet.legend());
          return 1;
        }
      }
    } else if (dset->DS()->Type() == DataSet::VECTOR) {
      // ----- Vector set --------------
      VarArray const& set_vars = *it;
      DataSet_Vector const& vecSet = static_cast<DataSet_Vector const&>( *(dset->DS()) );
      size_t start[2], count[2];
      start[1] = 0;
      count[0] = 1;
      count[1] = 3;
      for (unsigned int idx = 0; idx < vecSet.Size(); idx++) {
        start[0] = idx;
        if (NC::CheckErr(nc_put_vara(ncid_, set_vars[0].VID(), start, count, vecSet[idx].Dptr()))) {
          mprinterr("Error: Could not write vector %u xyz from '%s'\n", idx, vecSet.legend());
          return 1;
        }
      }
      if (vecSet.HasOrigins()) {
        for (unsigned int idx = 0; idx < vecSet.Size(); idx++) {
          start[0] = idx;
          if (NC::CheckErr(nc_put_vara(ncid_, set_vars[1].VID(), start, count, vecSet.OXYZ(idx).Dptr()))) {
            mprinterr("Error: Could not write vector %u origin from '%s'\n", idx, vecSet.legend());
            return 1;
          }
        }
      }
    } else if (dset->DS()->Type() == DataSet::XYMESH) {
      // ----- XY mesh set -------------
      VarArray const& set_vars = *it;
      // Write the X and Y variables
      size_t start[1];
      size_t count[1];
      start[0] = 0;
      count[0] = lengthDim.Size();

      DataSet_Mesh const& ds1d = static_cast<DataSet_Mesh const&>( *(dset->DS()) );
      std::vector<double> const& idxs = ds1d.MeshX();
      if (NC::CheckErr(nc_put_vara(ncid_, set_vars[1].VID(), start, count, &idxs[0]))) {
        mprinterr("Error: Could not write X variable from '%s'\n", ds1d.legend());
        return 1;
      }
      if (NC::CheckErr(nc_put_vara(ncid_, set_vars[0].VID(), start, count, ds1d.DvalPtr()))) {
        mprinterr("Error: Could not write variable '%s'\n", ds1d.legend());
        return 1;
      }
    } else {
      // ----- All other 1D sets -------
      size_t start[1];
      size_t count[1];
      start[0] = 0;
      count[0] = lengthDim.Size();
      DataSet_1D const& ds1d = static_cast<DataSet_1D const&>( *(dset->DS()) );
      if (NC::CheckErr(nc_put_vara(ncid_, it->front().VID(), start, count, ds1d.DvalPtr()))) {
        mprinterr("Error: Could not write variable '%s'\n", ds1d.legend());
        return 1;
      }
    }
  }
  return 0;
}

/** Write a 2D set. */
int DataIO_NetCDF::writeData_2D(DataSet const* ds) {
  // Define the dimension of the underlying array. Ensure name is unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  DataSet_2D const& set = static_cast<DataSet_2D const&>( *ds );
  int dimIdx = defineDim( "size", set.Size(), set.Meta().Legend() );
  if (dimIdx < 0) return 1;
  //mprintf("DEBUG: ncdim.Size()= %u\n", dimLen(dimIdx));
  // Choose type
  nc_type dtype;
  switch (set.Type()) {
    case DataSet::MATRIX_DBL : dtype = NC_DOUBLE ; break;
    case DataSet::MATRIX_FLT : dtype = NC_FLOAT ; break; 
    default:
      mprinterr("Internal Error: Unhandled DataSet type for 2D NetCDF variable.\n");
      return 1;
  }
  // Define the matrix variable
  NcVar matVar = defineVar(Dimensions_[dimIdx].DID(), dtype, ds->Meta().PrintName(), "matrix");
  if ( matVar.Empty() ) {
    mprinterr("Error: Could not define matrix variable for set '%s'\n", ds->legend());
    return 1;
  }
  // Add DataSet info
  if (AddDataSetInfo( ds, ncid_, matVar.VID() )) return 1;
  // Store rows and columns
  if (AddDataSetIntAtt( set.Ncols(), "ncols", ncid_, matVar.VID() )) return 1;
  if (AddDataSetIntAtt( set.Nrows(), "nrows", ncid_, matVar.VID() )) return 1;
  // Store the matrix kind
  std::string kind;
  switch (set.MatrixKind()) {
    case DataSet_2D::FULL : kind.assign("full"); break;
    case DataSet_2D::HALF : kind.assign("half"); break;
    case DataSet_2D::TRI  : kind.assign("tri"); break;
  }
  if (AddDataSetStringAtt(kind, "matrixkind", ncid_, matVar.VID())) return 1;
  // Define the diagonal vector/mass array if present
  NcVar vectVar, massVar;
  if (set.Type() == DataSet::MATRIX_DBL) {
    DataSet_MatrixDbl const& dmat = static_cast<DataSet_MatrixDbl const&>( set );
    if (dmat.Nsnapshots() > 0) {
      if (AddDataSetIntAtt( dmat.Nsnapshots(), "nsnapshots", ncid_, matVar.VID())) return 1;
    }
    if (!dmat.Vect().empty()) {
      int vDimIdx = defineDim("vectsize", dmat.Vect().size(), set.Meta().Legend() + " diagonal vector");
      if (vDimIdx < 0) return 1;
      vectVar = defineVar(Dimensions_[vDimIdx].DID(), NC_DOUBLE, ds->Meta().PrintName(), "vect", matVar.VID());
      if (vectVar.Empty()) { 
        mprinterr("Error: Could not define vect variable for matrix '%s'\n", set.legend());
        return 1;
      }
    }
    if (!dmat.Mass().empty()) {
      int mDimIdx = defineDim( "nmass", dmat.Mass().size(), set.Meta().Legend() + " mass" );
      if (mDimIdx < 0) return 1;
      massVar = defineVar(Dimensions_[mDimIdx].DID(), NC_DOUBLE, ds->Meta().PrintName(), "mass", matVar.VID());
      if (massVar.Empty()) {
        mprinterr("Error: Could not define mass variable for matrix '%s'\n", set.legend());
        return 1;
      }
    }
  }
  // END define variable
  if (EndDefineMode( ncid_ )) return 1;
  // Write the matrix 
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(dimIdx);
  //mprintf("DEBUG: start %zu count %zu\n", start[0], count[0]);
  if (NC::CheckErr(nc_put_vara(ncid_, matVar.VID(), start, count, set.MatrixPtr()))) {
    mprinterr("Error: Could not write matrix '%s'\n", set.legend());
    return 1;
  }
  if (!vectVar.Empty()) {
    DataSet_MatrixDbl const& dmat = static_cast<DataSet_MatrixDbl const&>( set );
    count[0] = dmat.Vect().size();
    if (NC::CheckErr(nc_put_vara(ncid_, vectVar.VID(), start, count, &(dmat.Vect()[0])))) {
      mprinterr("Error: Could not write vect variable for matrix '%s'\n", set.legend());
      return 1;
    }
  }
  if (!massVar.Empty()) {
    DataSet_MatrixDbl const& dmat = static_cast<DataSet_MatrixDbl const&>( set );
    count[0] = dmat.Mass().size();
    if (NC::CheckErr(nc_put_vara(ncid_, massVar.VID(), start, count, &(dmat.Mass()[0])))) {
      mprinterr("Error: Could not write mass variable for matrix '%s'\n", set.legend());
      return 1;
    }
  }

  return 0;
}

/** Write a 3D set. */
int DataIO_NetCDF::writeData_3D(DataSet const* ds) {
  // Define the dimension of the underlying array. Ensure name is unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  DataSet_3D const& set = static_cast<DataSet_3D const&>( *ds );
  int dimIdx = defineDim( "size", set.Size(), set.Meta().Legend() );
  if (dimIdx < 0) return 1;
  // Choose type
  nc_type dtype;
  switch (set.Type()) {
    case DataSet::GRID_DBL : dtype = NC_DOUBLE ; break;
    case DataSet::GRID_FLT : dtype = NC_FLOAT ; break; 
    default:
      mprinterr("Internal Error: Unhandled DataSet type for 3D NetCDF variable.\n");
      return 1;
  }
  // Define the grid variable
  NcVar gridVar = defineVar(Dimensions_[dimIdx].DID(), dtype, ds->Meta().PrintName(), "grid");
  if ( gridVar.Empty() ) {
    mprinterr("Error: Could not define grid variable for set '%s'\n", ds->legend());
    return 1;
  }
  // Add DataSet info to variable
  if (AddDataSetInfo( ds, ncid_, gridVar.VID() )) return 1;
  // Store NX, NY, and NZ
  if (AddDataSetIntAtt( set.NX(), "nx", ncid_, gridVar.VID() )) return 1;
  if (AddDataSetIntAtt( set.NY(), "ny", ncid_, gridVar.VID() )) return 1;
  if (AddDataSetIntAtt( set.NZ(), "nz", ncid_, gridVar.VID() )) return 1;
  // Store the origin
  if (NC::CheckErr(nc_put_att_double(ncid_, gridVar.VID(), "origin", NC_DOUBLE, 3,
                                     set.Bin().GridOrigin().Dptr() )))
  {
    mprinterr("Error: Could not put attribute 'origin'.\n");
    return 1;
  }
  // Store the unit cell vectors
  if (NC::CheckErr(nc_put_att_double(ncid_, gridVar.VID(), "ucell", NC_DOUBLE, 9,
                                     set.Bin().GridBox().UnitCell().Dptr() )))
  {
    mprinterr("Error: Could not put attribute 'ucell'.\n");
    return 1;
  }
  // END define variable
  if (EndDefineMode( ncid_ )) return 1;
  // Write the grid 
  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = dimLen(dimIdx);
  //mprintf("DEBUG: start %zu count %zu\n", start[0], count[0]);
  if (NC::CheckErr(nc_put_vara(ncid_, gridVar.VID(), start, count, set.GridPtr()))) {
    mprinterr("Error: Could not write grid '%s'\n", set.legend());
    return 1;
  }
  return 0;
}

/** Write pairwise cluster matrix to file. TODO fold into 2d matrix write? */
int DataIO_NetCDF::writeData_cluster_pwmatrix(DataSet const* ds) {
  // Define the dimensions of all arrays. Ensure names are unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  DataSet_PairwiseCache_MEM const& Mat = static_cast<DataSet_PairwiseCache_MEM const&>( *ds );

  // Actual matrix size
  int dimIdx = defineDim( "msize", Mat.Size(), Mat.Meta().Legend() + " matrix size" );
  if (dimIdx < 0) return 1;
  NcDim matsizeDim = Dimensions_[dimIdx];
  // Define variables TODO units
  // Matrix variable
  NcVar matrixVar = defineVar(matsizeDim.DID(), NC_FLOAT, Mat.Meta().PrintName(), "matrix");
  if (matrixVar.Empty()) {
    mprinterr("Error: Could not define cluster matrix variable.\n");
    return 1;
  }
  // Frames variable (if sieved) TODO always define?
  NcDim nrowsDim;
  NcVar framesVar;
  std::vector<int> indexVarIds;
  if (Mat.SieveVal() != 1) {
    // # matrix Rows (frames)
    dimIdx = defineDim( "n_rows", Mat.Nrows(), Mat.Meta().Legend() + " number of rows" );
    if (dimIdx < 0) return 1;
    nrowsDim = Dimensions_[dimIdx];
    // Actual frames clustered
    framesVar = defineVar(nrowsDim.DID(), NC_INT, Mat.Meta().PrintName(), "actual_frames");
    if (framesVar.Empty()) {
      mprinterr("Error: Could not define cluster frames variable.\n");
      return 1;
    }
    indexVarIds.assign(2, framesVar.VID());
  }
  // Add DataSet info to matrix variable
  if (AddDataSetInfo( ds, indexVarIds, ncid_, matrixVar.VID() )) return 1;
  // Attributes for matrix
  if (AddDataSetIntAtt( Mat.Nrows(), "ncols", ncid_, matrixVar.VID() )) return 1;
  if (AddDataSetIntAtt( Mat.Nrows(), "nrows", ncid_, matrixVar.VID() )) return 1;
  if (AddDataSetStringAtt("tri", "matrixkind", ncid_, matrixVar.VID())) return 1;
  if (AddDataSetIntAtt( Mat.FrameToIdx().size(), "n_original_frames", ncid_, matrixVar.VID() )) return 1;
  if (AddDataSetIntAtt( Mat.SieveVal(), "sieve", ncid_, matrixVar.VID() )) return 1;
  if (AddDataSetStringAtt(Mat.MetricDescrip(), "MetricDescription", ncid_, matrixVar.VID()))
    return 1;


  // END define variable
  if (EndDefineMode( ncid_ )) return 1;

  // Write variables
  size_t start[1], count[1];
  start[0] = 0;
  count[0] = matsizeDim.Size();
  if (NC::CheckErr(nc_put_vara(ncid_, matrixVar.VID(), start, count, Mat.Ptr()))) {
    mprinterr("Error: Writing cluster matrix variable.\n");
    return 1;
  }
  if (!framesVar.Empty()) {
    count[0] = nrowsDim.Size();
    Cpptraj::Cluster::Cframes actualFrames = Mat.PresentFrames();
    if (NC::CheckErr(nc_put_vara( ncid_, framesVar.VID(), start, count, &actualFrames[0] ))) {
      mprinterr("Error: Writing cluster frames variable.\n");
      return 1;
    }
  }
 
  return 0;
}

/** Write modes set to file. */
int DataIO_NetCDF::writeData_modes(DataSet const* ds) {
  // Define the dimensions of all arrays. Ensure names are unique by appending an index.
  if (EnterDefineMode(ncid_)) return 1;
  DataSet_Modes const& modes = static_cast<DataSet_Modes const&>( *ds );
  // N modes
  int dimIdx = defineDim( "nmodes", modes.Nmodes(), modes.Meta().Legend() + " eigenvalues" );
  if (dimIdx < 0) return 1;
  NcDim modesDim = Dimensions_[dimIdx];
  // N evec elements
  dimIdx = defineDim( "nevecelts", modes.Nmodes()*modes.VectorSize(), modes.Meta().Legend() + " eigenvectors" );
  if (dimIdx < 0) return 1;
  NcDim evecsDim = Dimensions_[dimIdx];
  // N avg coords
  dimIdx = defineDim( "ncoords", modes.NavgCrd(), modes.Meta().Legend() + " avg. coords" );
  if (dimIdx < 0) return 1;
  NcDim coordsDim = Dimensions_[dimIdx];
  // N Masses
  NcDim massDim;
  if (!modes.Mass().empty()) {
    dimIdx = defineDim( "nmass", modes.Mass().size(), modes.Meta().Legend() + " mass" );
    if (dimIdx == -1) return 1;
    massDim = Dimensions_[dimIdx];
  }
  // Define variables
  // Define the eigenvalues variable
  NcVar valuesVar = defineVar(modesDim.DID(), NC_DOUBLE, modes.Meta().PrintName(), "eigenvalues");
  if (valuesVar.Empty()) return 1;
  // Add DataSet info to variable
  if (AddDataSetInfo( ds, ncid_, valuesVar.VID() )) return 1;

  // Define the eigenvectors variable
  NcVar vectorsVar = defineVar(evecsDim.DID(), NC_DOUBLE, modes.Meta().PrintName(), "eigenvectors", valuesVar.VID());
  if (vectorsVar.Empty()) return 1;

  // Add the avg. coords variable
  NcVar coordsVar = defineVar(coordsDim.DID(), NC_DOUBLE, modes.Meta().PrintName(), "avgcoords", valuesVar.VID());
  if (coordsVar.Empty()) return 1;

  // Add the masses variable
  NcVar massVar;
  if (!massDim.Empty()) {
    massVar = defineVar(massDim.DID(), NC_DOUBLE, modes.Meta().PrintName(), "mass", valuesVar.VID());
    if (massVar.Empty()) return 1;
  }
  // END define variable
  if (EndDefineMode( ncid_ )) return 1;

  size_t start[1];
  size_t count[1];
  // Write the eigenvalues
  start[0] = 0;
  count[0] = modes.Nmodes();
  if (NC::CheckErr(nc_put_vara(ncid_, valuesVar.VID(), start, count, modes.EigenvaluePtr()))) {
    mprinterr("Error: Could not write eigenvalues from '%s'\n", modes.legend());
    return 1;
  }
  // Write the eigenvectors
  count[0] = modes.Nmodes()*modes.VectorSize();
  if (NC::CheckErr(nc_put_vara(ncid_, vectorsVar.VID(), start, count, modes.Eigenvectors()))) {
    mprinterr("Error: Could not write eigenvectors from '%s'\n", modes.legend());
    return 1;
  }
  // Write the avg coords
  count[0] = modes.NavgCrd();
  if (NC::CheckErr(nc_put_vara(ncid_, coordsVar.VID(), start, count, &(modes.AvgCrd()[0])))) {
    mprinterr("Error: Could not write avg. coords from '%s'\n", modes.legend());
    return 1;
  }
  // Write masses
  if (!massVar.Empty()) {
    count[0] = modes.Mass().size();
    if (NC::CheckErr(nc_put_vara(ncid_, massVar.VID(), start, count, &(modes.Mass()[0])))) {
      mprinterr("Error: Could not write mass from '%s'\n", modes.legend());
      return 1;
    }
  }

  return 0;
}
#endif /* BINTRAJ */

/// \return true if this is a valid 1D set to write
static inline bool isValid1dSetToWrite(DataSet const* ds) {
  return (ds->Group() == DataSet::SCALAR_1D ||
          ds->Type() == DataSet::STRING ||
          ds->Type() == DataSet::VECTOR
         );
}

// DataIO_NetCDF::WriteData()
int DataIO_NetCDF::WriteData(FileName const& fname, DataSetList const& dsl)
{
# ifdef BINTRAJ
  ncid_ = -1;
  Dimensions_.clear();
  // TODO check existing file
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid_ ) ))
    return 1;

  // Place incoming DataSets into a pool. As they are handled they will
  // be removed from the pool.
  SetPool setPool( dsl );

  // Check our incoming data sets. Try to find common dimensions.
  for (unsigned int idx = 0; idx < setPool.Nsets(); idx++)
  {
    if (setPool.IsUsed(idx)) continue;

    DataSet const* ds = setPool.Set( idx );
    if (debug_ > 0)
      mprintf("DEBUG: '%s'\n", ds->legend());
    if (ds->Type() == DataSet::MODES) {
      // ----- Modes -----------------------------
      if (writeData_modes(ds)) {
        mprinterr("Error: modes set write failed.\n");
        return 1;
      }
      setPool.MarkUsed( idx );
    } else if (ds->Type() == DataSet::PMATRIX_MEM) {
      // ----- Cluster pairwise matrix -----------
      if (writeData_cluster_pwmatrix(ds)) {
        mprinterr("Error: cluster pairwise matrix write failed.\n");
        return 1;
      }
      setPool.MarkUsed( idx );
    } else if (ds->Group() == DataSet::MATRIX_2D) {
      // ----- Matrix ----------------------------
      if (writeData_2D(ds)) {
        mprinterr("Error: matrix set write failed.\n");
        return 1;
      }
      setPool.MarkUsed( idx );
    } else if (ds->Group() == DataSet::GRID_3D) {
      // ----- Grid ------------------------------
      if (writeData_3D(ds)) {
        mprinterr("Error: grid set write failed.\n");
        return 1;
      }
    } else if (ds->Type() == DataSet::VECTOR_SCALAR) {
      // ----- Vector+Scalar ---------------------
      // TODO The VECTOR_SCALAR set has a dimension of 0 so it can default
      // to output with DataIO_Peaks, but this is a hack. Should be able to
      // be written with other 1D sets.
      if (writeData_vector_scalar(ds)) {
        mprinterr("Error: Vector/scalar set write failed.\n");
        return 1;
      }
    } else if (isValid1dSetToWrite(ds)) {
      // ----- 1D scalar -------------------------
      SetArray sets(1, Set(ds));
      setPool.MarkUsed( idx );
      // Group this set with others that share the same dimension and length.
      Dimension const& dim = ds->Dim(0);
      for (unsigned int jdx = idx + 1; jdx < setPool.Nsets(); jdx++)
      {
        DataSet const* ds2 = setPool.Set( jdx );
        if (isValid1dSetToWrite(ds2)) {
          if (ds->Size() == ds2->Size() && dim == ds2->Dim(0)) {
            sets.push_back( Set(ds2) );
            setPool.MarkUsed( jdx );
          }
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
    setPool.PrintUnused();
  }

  // Attributes
  if (EnterDefineMode(ncid_)) return 1;
  //if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"title",title.size(),title.c_str())) ) {
  //  mprinterr("Error: Writing title.\n");
  //  return 1;
  //}
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
  if (NC::PutConventions(ncid_, NC::NC_CPPTRAJDATA))
    return 1;
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"ConventionVersion",3,"1.0")) ) {
    mprinterr("Error: Writing NetCDF data conventions version.\n");
    return 1;
  }

  // End netcdf definitions
  if (EndDefineMode(ncid_)) {
    mprinterr("NetCDF data error on ending definitions.");
    return 1;
  }
  if (debug_ > 0)
    NC::Debug(ncid_);

  nc_close(ncid_); 
  return 0;
# else
  return 1;
# endif
}
