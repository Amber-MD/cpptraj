#include "Traj_H5.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Topology.h"
#ifdef BINTRAJ
# include <netcdf.h>
# include "NC_Routines.h"
#endif
//#ifdef HAS_HDF5
//# include <H5Cpp.h>
//  using namespace H5;
//#endif

/// CONSTRUCTOR
Traj_H5::Traj_H5()
//#ifdef HAS_HDF5
:
//file_(0)
  ncid_(-1),
  natom_(0)
//#endif
{}

/** DESTRUCTOR */
Traj_H5::~Traj_H5() {
  closeTraj();
//# ifdef HAS_HDF5
//  if (file_ != 0) delete file_;
//# endif
}

#ifdef HAS_HDF5
bool Traj_H5::HasConventions(int ncid) {
  std::string attrText = NC::GetAttrText(ncid, "conventions");
  if (attrText.empty())
    return false;
  if (attrText != "Pande") {
    mprinterr("Error: H5 file has conventions string that is not 'Pande'.\n");
    return false;
  }
  attrText = NC::GetAttrText(ncid, "conventionVersion");
  if ( attrText != "1.1")
    mprintf("Warning: H5 file has conventionVersion that is not 1.1 (%s)\n", attrText.c_str());
  return true;
}
#endif

/** Identify trajectory format. File should be setup for READ */
bool Traj_H5::ID_TrajFormat(CpptrajFile& fileIn) {
# ifdef HAS_HDF5
  int myNcid;
  if ( nc_open( fileIn.Filename().full(), NC_NOWRITE, &myNcid ) != NC_NOERR )
    return false;
  if (HasConventions(myNcid)) return true;
  nc_close( myNcid );
# else
  unsigned char buf[8];
  unsigned int nread = fileIn.Read(buf, 8);
  if (nread > 7 && buf[0] == 0x89 && buf[1] == 0x48 && buf[2] == 0x44 && buf[3] == 0x46 &&
                   buf[4] == 0x0d && buf[5] == 0x0a && buf[6] == 0x1a && buf[7] == 0x0a)
  {
    mprintf("Warning: File '%s' appears to be HDF5 but cpptraj was compiled without HDF5 support.\n", fname);
  }
# endif
  return false;
}

/** Print trajectory info to stdout. */
void Traj_H5::Info() {
  mprintf("is a MDtraj H5 (HDF5) trajectory");
}

/** Close file. */
void Traj_H5::closeTraj() {
# ifdef HAS_HDF5
  if (ncid_ == -1) return;
  bool err = NC::CheckErr( nc_close(ncid_) );
  if (err) {
    mprinterr("Error closing ncid %i\n", ncid_);
  }
  //if (ncdebug_ > 0 && !err)
  //  mprintf("Successfully closed ncid %i\n",ncid_);
  ncid_ = -1;
# endif
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_H5::openTrajin() {

  return 0;
}

/** Read help */
void Traj_H5::ReadHelp() {

}

/** Process read arguments. */
int Traj_H5::processReadArgs(ArgList& argIn) {

  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_H5::setupTrajin(FileName const& fname, Topology* trajParm)
{
# ifdef HAS_HDF5
/*  if (file_ != 0)
    delete file_;
  file_ = new H5File( fname.full(), H5F_ACC_RDONLY );
  if (file_ == 0) {
    mprinterr("Error: Could not allocate input trajectory file.\n");
    return TRAJIN_ERR;
  }
  return 0;*/
  if (ncid_ != -1) closeTraj();
  if ( NC::CheckErr( nc_open( fname.full(), NC_NOWRITE, &ncid_ ) ) )
    return TRAJIN_ERR;
  NC::Debug(ncid_);
  // get variables/dims present
  // ndimsp:  Pointer to location for returned number of dimensions defined for 
  //          this netCDF dataset.
  // nvarsp:  Pointer to location for returned number of variables defined for 
  //          this netCDF dataset.
  // ngattsp: Pointer to location for returned number of global attributes 
  //          defined for this netCDF dataset.
  // unlimdimidp: 
  //  Pointer to location for returned ID of the unlimited dimension, if 
  //  there is one for this netCDF dataset. If no unlimited length 
  //  dimension has been defined, -1 is returned
  int ndimsp, nvarsp, ngattsp, unlimdimidp;
  char varname[NC_MAX_NAME+1];
  if (NC::CheckErr( nc_inq(ncid_, &ndimsp, &nvarsp, &ngattsp, &unlimdimidp) ) )
    return TRAJIN_ERR;
  mprintf("DEBUG: Unlimited dimid is %i\n", unlimdimidp);
  if (unlimdimidp < 0) {
    mprinterr("Error: No unlimited dimension present in H5 file.\n");
    return TRAJIN_ERR;
  }
  // Check for recognized variables.
  int nframes = 0;
  for (int ivar = 0; ivar < nvarsp; ivar++) {
    if (NC::CheckErr( nc_inq_varname(ncid_, ivar, varname) ) )
      return TRAJIN_ERR;
    std::string varstr( varname );
    if (varstr == "coordinates") {
      // Get dimensions
      int ndims = 0;
      if (NC::CheckErr( nc_inq_varndims(ncid_, ivar, &ndims) ) )
        return TRAJIN_ERR;
      if (ndims != 3) {
        mprinterr("Error: Expected 3 dims for 'coordinates', got %i\n", ndims);
        return TRAJIN_ERR;
      }
      int coord_dims[3];
      if (NC::CheckErr( nc_inq_vardimid(ncid_, ivar, coord_dims) ) )
        return TRAJIN_ERR;
      mprintf("DEBUG: Coord dims: %i %i %i\n", coord_dims[0], coord_dims[1], coord_dims[2]);
      // Check the dimensions. One should be frames (unlimited), one should be
      // atoms, and the last should be spatial (XYZ)
      bool dim_used[3];
      for (int i = 0; i < 3; i++) dim_used[i] = false;
      bool has_unlimited = false;
      natom_ = 0;
      for (int nd = 0; nd < ndims; nd++) {
        size_t dimsize = 0;
        if (NC::CheckErr( nc_inq_dimlen(ncid_, coord_dims[nd], &dimsize)))
          return TRAJIN_ERR;
        mprintf("DEBUG: Dim %i size %zu\n", coord_dims[nd], dimsize);
        if (coord_dims[nd] == unlimdimidp) {
          has_unlimited = true;
          dim_used[nd] = true;
          nframes = (int)dimsize;
        } else if ( dimsize == 3 ) {
          dim_used[nd] = true;
          mprintf("DEBUG: Spatial dim is %i\n", coord_dims[nd]);
        } else {
          // Should be natoms
          dim_used[nd] = true;
          natom_ = (int)dimsize;
          mprintf("DEBUG: Atom dim is %i\n", coord_dims[nd]);
        }
      }
      // Check # atoms
      if ( natom_ != trajParm->Natom() ) {
        mprinterr("Error: Atom mismatch between topology (%i) and trajectory (%i).\n",
                   trajParm->Natom(), natom_);
        return TRAJIN_ERR;
      }
      if (!has_unlimited) {
        mprinterr("Error: No unlimited dimension for 'coordinates'.\n");
        return TRAJIN_ERR;
      }
      // Expect 1 to be spatial, 2 to be natom
      // TODO is this always the case for H5?
      
      
    } else {
      mprintf("Warning: Unhandled variable present: %s\n", varname);
    }
  }

  return 0;
# else
  return TRAJIN_ERR;
# endif
}

/** Read specified trajectory frame. */
int Traj_H5::readFrame(int set, Frame& frameIn) {

  return 0;
}

/** Read velocities from specified frame. */
int Traj_H5::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_H5::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_H5::WriteHelp() {

}

/** Process write arguments. */
int Traj_H5::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_H5::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_H5::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_H5::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_H5::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_H5::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_H5::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_H5::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_H5::parallelCloseTraj() {

}
#endif
