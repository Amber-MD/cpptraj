#ifdef BINTRAJ
#include <netcdf.h>
#include "DataSet_integer_disk.h"
#include "CpptrajStdio.h"
#include "NC_Routines.h"
#include "File_TempName.h"

/// CONSTRUCTOR
DataSet_integer_disk::DataSet_integer_disk() :
  ncid_(-1),
  framevid_(-1),
  nvals_(0)
{
  start_[0] = 0;
  count_[0] = 0;
  // NOTE: Since tmpnam() use is considered dangerous, use the
  //       DataSetDiskCache name in its place.
  tfname_ = File::GenTempName();
  if (tfname_.empty()) {
    mprinterr("Internal Error: Could not get temporary file name of dist integer data set.\n");
    return;
  }
  if (NC::CheckErr( nc_create( tfname_.full(), NC_64BIT_OFFSET, &ncid_ ) )) {
    mprinterr("Internal Error: Could not disk cache integer data set.\n");
  }
  mprintf("DEBUG: Integer data set being cached in file: %s\n", tfname_.full());
  int frameDID;
  if (NC::CheckErr(nc_def_dim(ncid_, "frame", NC_UNLIMITED, &frameDID))) {
    mprinterr("Internal Error: Could not define frame dimension for disk integer data set.\n");
  }
  int dimensionID[NC_MAX_VAR_DIMS];
  dimensionID[0] = frameDID;
  if (NC::CheckErr(nc_def_var(ncid_, "values", NC_INT, 1, dimensionID, &framevid_))) {
    mprinterr("Internal Error: Could not define frame variable for disk integer data set.\n");
  }
  // End definitions
  if (NC::CheckErr(nc_enddef(ncid_))) {
    mprinterr("Internal Error: Ending definitions for disk integer data set.");
  }
}

/// DESTRUCTOR
DataSet_integer_disk::~DataSet_integer_disk() {
  if (ncid_ != -1) nc_close( ncid_ );
  if (!tfname_.empty()) File::FreeTempName( tfname_ );
} 

#ifdef MPI
// TODO
int DataSet_integer_disk::Sync(size_t total, std::vector<int> const& rank_frames,
                          Parallel::Comm const& commIn)
{
  if (commIn.Size() > 1) {
    mprinterr("Internal Error: Data set sync for integer disk-cached data sets is not yet supported.\n");
    return 1;
  }
  return 0;
}
#endif

// DataSet_integer_disk::Info()
void DataSet_integer_disk::Info() const {
  mprintf(" (cached)");
}

// DataSet_integer_disk::Allocate()
int DataSet_integer_disk::Allocate(SizeArray const& sizeIn) {
  return 0;
}

//  DataSet_integer_disk::getVal()
int DataSet_integer_disk::getVal(size_t idx) const {
  size_t start[1], count[1];
  start[0] = idx;
  count[0] = 1;
  int val;
  nc_get_vara_int(ncid_, framevid_, start, count, &val);
  return val;
}

//  DataSet_integer_disk::setVal()
void DataSet_integer_disk::setVal(size_t idx, int val) {
  start_[0] = idx;
  count_[0] = 1;
  nc_put_vara_int(ncid_, framevid_, start_, count_, &val);
  //nc_sync(ncid_);
}

// DataSet_integer_disk::Add()
void DataSet_integer_disk::Add(size_t frame, const void* vIn) {
  if (frame > nvals_)
    Resize(frame);
  start_[0] = nvals_;
  count_[0] = 1;
  nc_put_vara_int(ncid_, framevid_, start_, count_, (const int*)vIn);
  //nc_sync(ncid_);
  nvals_++;
}

// DataSet_integer_disk::WriteBuffer()
void DataSet_integer_disk::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= nvals_)
    cbuffer.Printf(format_.fmt(), 0);
  else
    cbuffer.Printf(format_.fmt(), getVal(pIn[0]));
}

// DataSet_integer_disk::Append()
int DataSet_integer_disk::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != SCALAR_1D) return 1;
  DataSet_1D const& ds = static_cast<DataSet_1D const&>( *dsIn );
  for (unsigned int i = 0; i != ds.Size(); i++) {
    int val =(int)ds.Dval(i);
    Add(i, &val);
  }
  return 0;
}

// DataSet_integer_disk::Dval()
double DataSet_integer_disk::Dval(size_t idx) const {
  return (double)getVal(idx);
}

//  DataSet_integer_disk::VoidPtr()
const void* DataSet_integer_disk::VoidPtr(size_t offset) const {
  mprinterr("Internal Error: VoidPtr() not implemented for DataSet_integer_disk.\n");
  return 0;
}

// DataSet_integer_disk::SetElement()
void DataSet_integer_disk::SetElement(size_t idx, int val) {
  setVal(idx, val);
}

// DataSet_integer_disk::operator[]()
int DataSet_integer_disk::operator[](size_t idx) const {
  return getVal(idx);
}

// FIXME redundant
void DataSet_integer_disk::AddElement(int val) {
  Add(0, &val);
}

// DataSet_integer_disk::Resize()
void DataSet_integer_disk::Resize(size_t sizeIn) {
  if (sizeIn < nvals_)
    nvals_ = sizeIn;
  else if (sizeIn > nvals_) {
    const int zero = 0;
    for (size_t idx = nvals_; idx != sizeIn; idx++)
      setVal(idx, zero);
    nvals_ = sizeIn;
  }
}

// DataSet_integer_disk::Assign()
void DataSet_integer_disk::Assign(size_t sizeIn, int val) {
  for (size_t idx = 0; idx != sizeIn; idx++)
    setVal(idx, val);
  nvals_ = sizeIn;
}

// DataSet_integer_disk::AddVal()
void DataSet_integer_disk::AddVal(size_t idx, int val) {
  if (idx < nvals_) {
    int newval = getVal(idx) + val;
    setVal(idx, newval);
  } else
    Add(idx, &val);
}
#endif /* BINTRAJ */
