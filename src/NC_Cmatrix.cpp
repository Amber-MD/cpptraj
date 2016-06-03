#ifdef BINTRAJ
# include <netcdf.h>
# include "NC_Routines.h"
#endif
#include "NC_Cmatrix.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
NC_Cmatrix::NC_Cmatrix() :
  ncid_(-1),
  n_original_frames_DID_(-1),
  n_rows_DID_(-1),
  msize_DID_(-1),
  cmatrix_VID_(-1),
  actualFrames_VID_(-1),
  nRows_(0),
  mSize_(0)
{}

NC_Cmatrix::~NC_Cmatrix() {
  CloseCmatrix();
}

// NC_Cmatrix::ID_Cmatrix()
bool NC_Cmatrix::ID_Cmatrix(FileName const& fname) {
# ifdef BINTRAJ
  int NCID;
  if ( nc_open( fname.full(), NC_NOWRITE, &NCID ) != NC_NOERR )
    return false;
  std::string attrText = NC::GetAttrText(NCID, "Conventions");
  nc_close( NCID );
  return (attrText == "CPPTRAJ_CMATRIX");
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
int NC_Cmatrix::OpenCmatrixRead(FileName const& fname) {
  if (ncid_ != -1) CloseCmatrix();
  if (fname.empty()) return 1;
  if (NC::CheckErr( nc_open( fname.full(), NC_NOWRITE, &ncid_ ) ))
    return 1;
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

// NC_Cmatrix::OpenCmatrixWrite()
int NC_Cmatrix::OpenCmatrixWrite(FileName const& fname, unsigned int nFrames, unsigned int nRowsIn,
                                 int sieve)
{
  mprinterr("DEBUG: Cmatrix file '%s', nFrames %u, nRows %u, sieve %i\n",
          fname.full(), nFrames, nRowsIn, sieve);
  if (fname.empty()) return 1;
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid_ ) ))
    return 1;

  nRows_ = nRowsIn;
  if (nRows_ < 1) {
    mprinterr("Internal Error: Trying to create empty cluster matrix file.\n");
    return 1;
  }
  // Define dimensions
  if (NC::CheckErr( nc_def_dim( ncid_, NC_CMATRIX_NFRAMES, nFrames, &n_original_frames_DID_ ) ))
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
  //if (NC::CheckErr(nc_set_fill(ncid_, NC_NOFILL, dimensionID))) {
  //  mprinterr("Error: NetCDF setting fill value.\n");
  //  return 1;
  //}

  // End netcdf definitions
  if (NC::CheckErr(nc_enddef(ncid_))) return 1;

  // Write sieve value
  if (NC::CheckErr(nc_put_var_int(ncid_, sieveVID, &sieve))) return 1;

  return 0;
}

// NC_Cmatrix::WriteFramesArray()
int NC_Cmatrix::WriteFramesArray(std::vector<int> const& actualFrames) {
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
int NC_Cmatrix::WriteCmatrixElement(unsigned int xIn, unsigned int yIn, double dval)
{
  float fval = (float)dval;
  size_t index[1];
  long int idx = CalcIndex(xIn, yIn);
  if (idx < 0L) return 1;
  if (NC::CheckErr( nc_put_var1_float(ncid_, cmatrix_VID_, index, &fval) ))
    return 1;
  return 0;
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
int OpenCmatrixRead(FileName const&) { return 1; }
double GetCmatrixElement(unsigned int, unsigned int) { return 0.0; }
int OpenCmatrixWrite(FileName const&, int, int, int) { return 1; }
int WriteFramesArray(std::vector<int> const&) { return 1; }
int WriteCmatrixElement(unsigned int, unsigned int, double) { return 1; }
void CloseCmatrix() {}
#endif
