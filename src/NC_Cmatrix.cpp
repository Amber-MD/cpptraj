#ifdef BINTRAJ
# include <netcdf.h>
# include "NC_Routines.h"
#endif
#include "NC_Cmatrix.h"
#include "CpptrajStdio.h"

#ifdef BINTRAJ
/// CONSTRUCTOR
NC_Cmatrix::NC_Cmatrix() :
  ncid_(-1),
  n_original_frames_DID_(-1),
  n_rows_DID_(-1),
  msize_DID_(-1),
  cmatrix_VID_(-1),
  actualFrames_VID_(-1),
  nFrames_(0),
  nRows_(0),
  mSize_(0),
  mode_(READ)
{}

NC_Cmatrix::~NC_Cmatrix() {
  CloseCmatrix();
}

bool NC_Cmatrix::IsCpptrajCmatrix(int NCID) {
  return (NC::GetAttrText(NCID, "Conventions") == "CPPTRAJ_CMATRIX");
}
#endif

// NC_Cmatrix::ID_Cmatrix()
bool NC_Cmatrix::ID_Cmatrix(FileName const& fname) {
# ifdef BINTRAJ
  int NCID;
  if ( nc_open( fname.full(), NC_NOWRITE, &NCID ) != NC_NOERR )
    return false;
  bool isCmatrix = IsCpptrajCmatrix(NCID);
  nc_close( NCID );
  return isCmatrix;
# else
  mprinterr("Error: Compiled without NetCDF support. Recompile with -DBINTRAJ.\n");
  return false;
# endif
}

#ifdef BINTRAJ
// DEFINES
#define NC_CMATRIX_NFRAMES "n_original_frames"
#define NC_CMATRIX_NROWS "n_rows"
#define NC_CMATRIX_MSIZE "msize"
#define NC_CMATRIX_SIEVE "sieve"
#define NC_CMATRIX_MATRIX "matrix"
#define NC_CMATRIX_FRAMES "actual_frames"

// NC_Cmatrix::OpenCmatrixRead()
int NC_Cmatrix::OpenCmatrixRead(FileName const& fname, int& sieve) {
  if (ncid_ != -1) CloseCmatrix();
  if (fname.empty()) return 1;
  if (NC::CheckErr( nc_open( fname.full(), NC_NOWRITE, &ncid_ ) ))
    return 1;
  if (!IsCpptrajCmatrix(ncid_)) {
    mprinterr("Error: File '%s' is not cpptraj cluster matrix.\n", fname.full());
    return 1;
  }
  mode_ = READ;

  // Attributes
  std::string version = NC::GetAttrText(ncid_, "Version");
  if (version != "1.0")
    mprintf("Warning: NetCDF cluster matrix file is version '%s'; expected '1.0'\n",
            version.c_str());

  // Dimensions
  n_original_frames_DID_ = NC::GetDimInfo(ncid_, NC_CMATRIX_NFRAMES, nFrames_);
  if (n_original_frames_DID_ == -1) {
    mprinterr("Error: Could not get frames dimension.\n");
    return 1;
  }
  n_rows_DID_ = NC::GetDimInfo(ncid_, NC_CMATRIX_NROWS, nRows_);
  if (n_rows_DID_ == -1) {
    mprinterr("Error: Could not get rows dimension.\n");
    return 1;
  }
  msize_DID_ = NC::GetDimInfo(ncid_, NC_CMATRIX_MSIZE, mSize_);
  if (msize_DID_ == -1) {
    mprinterr("Error: Could not get matrix size dimension.\n");
    return 1;
  }
  // Variables
  // Sieve
  int sieveVID;
  if (NC::CheckErr(nc_inq_varid(ncid_, NC_CMATRIX_SIEVE, &sieveVID))) {
    mprinterr("Error: Could not get sieve variable id.\n");
    return 1;
  }
  if (NC::CheckErr(nc_get_var_int(ncid_, sieveVID, &sieve))) return 1;
  // Matrix
  if (NC::CheckErr(nc_inq_varid(ncid_, NC_CMATRIX_MATRIX, &cmatrix_VID_))) {
    mprinterr("Error: Could not get matrix variable id.\n");
    return 1;
  }
  // Frames; only allowed to be -1 if sieve is 1
  if (nc_inq_varid(ncid_, NC_CMATRIX_FRAMES, &actualFrames_VID_) != NC_NOERR) {
    if (sieve != 1) {
      mprinterr("Error: Cluster matrix has sieve but no frames variable id.\n");
      return 1;
    }
    actualFrames_VID_ = -1;
  }

  return 0;
}

// NC_Cmatrix::CalcIndex()
// TODO Consolidate code in Matrix.h?
long int NC_Cmatrix::CalcIndex(unsigned int xIn, unsigned int yIn) const {
  // Calculate upper-triangle matrix index
  unsigned int i, j;
  if (yIn > xIn) {
    i = xIn;
    j = yIn;
  } else if (xIn > yIn) {
    i = yIn;
    j = xIn;
  } else { // iIn == jIn, triangle matrix diagonal is not valid
    mprinterr("Error: Invalid attempt to access diagonal from cluster matrix (%i,%i)\n", xIn, yIn);
    return -1L;
  }
  unsigned int i1 = i + 1;
  return (long int)( ( (nRows_ * i) - ((i1 * i) / 2) ) + j - i1 );
}

// NC_Cmatrix::GetCmatrixElement()
double NC_Cmatrix::GetCmatrixElement(unsigned int xIn, unsigned int yIn) const {
  float fval;
  size_t index[1];
  long int idx = CalcIndex(xIn, yIn);
  if (idx < 0L) return 0.0;
  index[0] = (size_t)idx;
  if (NC::CheckErr( nc_get_var1_float(ncid_, cmatrix_VID_, index, &fval) ))
    return 0.0;
  return (double)fval;
}

// NC_Cmatrix::GetCmatrixElement()
double NC_Cmatrix::GetCmatrixElement(unsigned int idxIn) const {
  float fval;
  size_t index[1] = { idxIn };
  if (NC::CheckErr( nc_get_var1_float(ncid_, cmatrix_VID_, index, &fval) ))
    return 0.0;
  return (double)fval;
}

// NC_Cmatrix::GetSieveStatus()
std::vector<char> NC_Cmatrix::GetSieveStatus() const {
  if (nFrames_ < 1) return std::vector<char>();
  if (actualFrames_VID_ == -1)
    // No frames array. All frames present, none sieved.
    return std::vector<char>(nFrames_, 'F');
  else {
    // Get the frames array
    std::vector<int> actualFrames(nRows_);
    size_t start[1] = { 0      };
    size_t count[1] = { nRows_ };
    if (NC::CheckErr(nc_get_vara_int( ncid_, actualFrames_VID_, start, count, &actualFrames[0] )))
      return std::vector<char>();
    // All frames execpt those in actualFrames are sieved out.
    std::vector<char> sieveStatus(nFrames_, 'T');
    for (std::vector<int>::const_iterator it = actualFrames.begin();
                                          it != actualFrames.end(); ++it)
      sieveStatus[ *it ] = 'F';
    return sieveStatus;
  }
}

// NC_Cmatrix::GetCmatrix()
int NC_Cmatrix::GetCmatrix(float* ptr) const {
  if (cmatrix_VID_ == -1) return 1;
  size_t start[1] = { 0      };
  size_t count[1] = { mSize_ };
  return NC::CheckErr(nc_get_vara_float( ncid_, cmatrix_VID_, start, count, ptr ));
}

// NC_Cmatrix::CreateCmatrix()
int NC_Cmatrix::CreateCmatrix(FileName const& fname, unsigned int nFramesIn, unsigned int nRowsIn,
                              int sieve)
{
  //mprinterr("DEBUG: Cmatrix file '%s', nFrames %u, nRows %u, sieve %i\n",
  //        fname.full(), nFrames, nRowsIn, sieve);
  if (fname.empty()) return 1;
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid_ ) ))
    return 1;
  nFrames_ = nFramesIn;
  nRows_ = nRowsIn;
  if (nRows_ < 1) {
    mprinterr("Internal Error: Trying to create empty cluster matrix file.\n");
    return 1;
  }
  mode_ = WRITE;
  // Define dimensions
  if (NC::CheckErr( nc_def_dim( ncid_, NC_CMATRIX_NFRAMES, nFrames_, &n_original_frames_DID_ ) ))
    return 1;
  if (NC::CheckErr( nc_def_dim( ncid_, NC_CMATRIX_NROWS, nRows_, &n_rows_DID_ ) ))
    return 1;
  mSize_ = (nRows_ * (nRows_-1)) / 2;
  if (NC::CheckErr( nc_def_dim( ncid_, NC_CMATRIX_MSIZE, mSize_, &msize_DID_ ) ))
    return 1;

  // Define variables // TODO units
  int dimensionID[1];
  // Sieve
  int sieveVID;
  if (NC::CheckErr(nc_def_var( ncid_, NC_CMATRIX_SIEVE, NC_INT, 0, dimensionID, &sieveVID ) )) {
    mprinterr("Error: Defining sieve variable.\n");
    return 1;
  }
  // Matrix
  dimensionID[0] = msize_DID_;
  if (NC::CheckErr(nc_def_var( ncid_, NC_CMATRIX_MATRIX, NC_FLOAT, 1, dimensionID,
                               &cmatrix_VID_ ) ))
  {
    mprinterr("Error: Defining matrix variable.\n");
    return 1;
  }
  // Frames (if sieved)
  if (sieve != 1) {
    dimensionID[0] = n_rows_DID_;
    if (NC::CheckErr(nc_def_var( ncid_, NC_CMATRIX_FRAMES, NC_INT, 1, dimensionID,
                                 &actualFrames_VID_ ) ))
    {
      mprinterr("Error: Defining actual frames variable.\n");
      return 1;
    }
  } else
    actualFrames_VID_ = -1;

  // Attributes
  if (NC::CheckErr(nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", 15, "CPPTRAJ_CMATRIX")))
    return 1;
  if (NC::CheckErr(nc_put_att_text(ncid_, NC_GLOBAL, "Version", 3, "1.0")))
    return 1;

  // Set fill mode
  if (NC::CheckErr(nc_set_fill(ncid_, NC_NOFILL, dimensionID))) {
    mprinterr("Error: NetCDF setting fill value.\n");
    return 1;
  }

  // End netcdf definitions
  if (NC::CheckErr(nc_enddef(ncid_))) return 1;

  // Write sieve value
  if (NC::CheckErr(nc_put_var_int(ncid_, sieveVID, &sieve))) return 1;

  return 0;
}

// NC_Cmatrix::Sync()
void NC_Cmatrix::Sync() const {
  if (ncid_ != -1)
    NC::CheckErr( nc_sync(ncid_) );
}

int NC_Cmatrix::ReopenSharedWrite(FileName const& fname) {
  if (ncid_ == -1) return 1;
  // Close and re-open shared
  nc_close( ncid_ );
  if (NC::CheckErr(nc_open(fname.full(), NC_WRITE|NC_SHARE, &ncid_))) return 1;
  return 0;
}

// NC_Cmatrix::WriteFramesArray()
int NC_Cmatrix::WriteFramesArray(std::vector<int> const& actualFrames) const {
  if (ncid_ == -1) return 1; // Sanity check
  if (actualFrames_VID_ == -1) {
    mprinterr("Error: No cluster frames variable ID defined.\n");
    return 1;
  }
  if (actualFrames.size() != nRows_) {
    mprinterr("Error: Frames array is %zu elements but expected %u\n",
              actualFrames.size(), nRows_);
    return 1;
  }
  size_t start[1] = { 0      };
  size_t count[1] = { nRows_ };
  if (NC::CheckErr(nc_put_vara_int( ncid_, actualFrames_VID_, start, count, &actualFrames[0] )))
    return 1;
  return 0;
}

// NC_Cmatrix::WriteCmatrixElement()
int NC_Cmatrix::WriteCmatrixElement(unsigned int xIn, unsigned int yIn, double dval) const
{
  int err = 0;
# ifdef _OPENMP
  // Since NetCDF files are not thread safe this must be protected.
# pragma omp critical(writecmatrixelement)
  {
# endif
  float fval = (float)dval;
  size_t index[1];
  long int idx = CalcIndex(xIn, yIn);
  //mprintf("DEBUG: Index calc x=%u y=%u i=%li\n", xIn, yIn, idx);
  if (idx < 0L)
    err = 1;
  else {
    index[0] = (size_t)idx;
    if (NC::CheckErr( nc_put_var1_float(ncid_, cmatrix_VID_, index, &fval) ))
      err = 1;
  }
# ifdef _OPENMP
  } // END pragma omp critical
# endif
  return err;
}

// NC_Cmatrix::WriteCmatrix()
int NC_Cmatrix::WriteCmatrix(const float* ptr) const {
  if (cmatrix_VID_ == -1) return 1;
  size_t start[1] = { 0      };
  size_t count[1] = { mSize_ };
  return NC::CheckErr(nc_put_vara_float( ncid_, cmatrix_VID_, start, count, ptr ));
}

// NC_Cmatrix::CloseCmatrix()
void NC_Cmatrix::CloseCmatrix() {
  if (ncid_ != -1) {
    nc_close( ncid_ );
    ncid_ = -1;
    // TODO: Set all IDs to -1 too?
  }
}
#else
NC_Cmatrix::NC_Cmatrix() {}
NC_Cmatrix::~NC_Cmatrix() {}
int NC_Cmatrix::OpenCmatrixRead(FileName const&, int&) { return 1; }
double NC_Cmatrix::GetCmatrixElement(unsigned int, unsigned int) const { return 0.0; }
double NC_Cmatrix::GetCmatrixElement(unsigned int) const { return 0.0; }
int NC_Cmatrix::CreateCmatrix(FileName const&, unsigned int, unsigned int, int) {
  mprinterr("Error: Cpptraj was compiled without NetCDF. Cannot create NetCDF matrix file.\n");
  return 1;
}
int NC_Cmatrix::WriteFramesArray(std::vector<int> const&) { return 1; }
int NC_Cmatrix::WriteCmatrixElement(unsigned int, unsigned int, double) { return 1; }
void NC_Cmatrix::CloseCmatrix() {}
void NC_Cmatrix::Sync() {}
int NC_Cmatrix::ReopenSharedWrite(FileName const&) { return 1; }
#endif
