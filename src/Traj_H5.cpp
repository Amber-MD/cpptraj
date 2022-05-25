#include "Traj_H5.h"
#include "Constants.h"
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
  natom_(0),
  coordVID_(-1),
  cellLengthVID_(-1),
  cellAngleVID_(-1),
  timeVID_(-1),
  convert_h5_to_cpptraj_box_(0),
  convert_h5_to_cpptraj_coord_(0)
//#endif
{
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 0;
  count_[1] = 0;
  count_[2] = 0;
}

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

# ifdef HAS_HDF5
static inline int setLengthFac(double& fac, std::string const& Units, const char* desc)
{
  if (Units == "nanometers") {
    fac = Constants::NM_TO_ANG;
  } else if (Units == "angstroms") {
    fac = 1.0;
  } else {
    mprinterr("Error: %s has unrecognized units (%s)\n", desc, Units.c_str());
    return 1;
  }
  return 0;
}

/** Set up the coordinates variable ID, number of atoms, and number of frames. */
int Traj_H5::setupCoordVID(int& frameDID, int& atomDID, int& spatialDID, int& nframes)
{
  // Get the 'coordinates' variable ID
  if (NC::CheckErr(nc_inq_varid(ncid_, "coordinates", &coordVID_)))
    return 1;
  mprintf("DEBUG: Coordinates VID is %i\n", coordVID_);
  // Set conversion factor for coords
  std::string lengthUnits = NC::GetAttrText(ncid_, coordVID_, "units");
  if (setLengthFac(convert_h5_to_cpptraj_coord_, lengthUnits, "Coordinates"))
    return 1;
  // Dimensions
  frameDID = -1;
  atomDID = -1;
  spatialDID = -1;
  natom_ = 0;
  nframes = 0;
  // Need to get the unlimited dimension ID, which should be the frame dim.
  if (NC::CheckErr( nc_inq_unlimdim(ncid_, &frameDID) ) )
    return 1;
  if (frameDID < 0) {
    mprinterr("Error: No unlimited (frame) dimension present in H5 file.\n");
    return 1;
  }
  // Get dimensions for coordinates
  int ndims = 0;
  if (NC::CheckErr( nc_inq_varndims(ncid_, coordVID_, &ndims) ) )
    return 1;
  if (ndims != 3) {
    mprinterr("Error: Expected 3 dims for 'coordinates', got %i\n", ndims);
    return 1;
  }
  int coord_dims[3];
  if (NC::CheckErr( nc_inq_vardimid(ncid_, coordVID_, coord_dims) ) )
    return 1;
  mprintf("DEBUG: Coord dims: %i %i %i\n", coord_dims[0], coord_dims[1], coord_dims[2]);
  // Check the dimensions. One should be frames (unlimited), one should be
  // atoms, and the last should be spatial (XYZ)

  bool dim_used[3];
  for (int i = 0; i < 3; i++) dim_used[i] = false;
  bool has_unlimited = false;
  for (int nd = 0; nd < ndims; nd++) {
    size_t dimsize = 0;
    if (NC::CheckErr( nc_inq_dimlen(ncid_, coord_dims[nd], &dimsize)))
      return 1;
    mprintf("DEBUG: Dim %i size %zu\n", coord_dims[nd], dimsize);
    if (coord_dims[nd] == frameDID) {
      has_unlimited = true;
      dim_used[nd] = true;
      nframes = (int)dimsize;
    } else if ( dimsize == 3 ) {
      dim_used[nd] = true;
      spatialDID = coord_dims[nd];
    } else {
      // Should be natoms
      dim_used[nd] = true;
      natom_ = (int)dimsize;
      atomDID = coord_dims[nd];
    }
  }

  // Sanity checks
  for (int nd = 0; nd < ndims; nd++) {
    if (!dim_used[nd]) {
      mprinterr("Error: Dimension %i remains unused for 'coordinates'.\n", coord_dims[nd]);
      return 1;
    }
  }
  if (!has_unlimited) {
    mprinterr("Error: No unlimited dimension for 'coordinates'.\n");
    return 1;
  }
  // Expect dim order to be frame, atom, spatial
  // TODO is this always the case for H5?
  if (coord_dims[0] != frameDID) {
    mprinterr("Error: Frame dimension is not first.\n");
    return 1;
  }
  if (coord_dims[1] != atomDID) {
    mprinterr("Error: Atom dimension is not second.\n");
    return 1;
  }
  if (coord_dims[2] != spatialDID) {
    mprinterr("Error: Spatial dimension is not third.\n");
    return 1;
  }
  return 0;
}

/** Set up box variable IDs.
  * \return 0 on success, 1 on error, -1 for no box coords.
  */
int Traj_H5::setupBoxVIDs(Box& ncbox, int frameDID, int spatialDID) {
  ncbox.SetNoBox();
  // Get the 'cell_lengths' variable ID
  int err = nc_inq_varid(ncid_, "cell_lengths", &cellLengthVID_);
  if (err != NC_NOERR) return -1;
  // Get the 'cell_angles' variable ID 
  if (NC::CheckErr(nc_inq_varid(ncid_, "cell_angles", &cellAngleVID_)))
    return 1;
  mprintf("DEBUG: Cell length vid= %i, cell angle vid= %i\n",
          cellLengthVID_, cellAngleVID_);
  // Ensure angles are in degrees
  std::string angleUnits = NC::GetAttrText(ncid_, cellAngleVID_, "units");
  if (angleUnits != "degrees") {
    mprinterr("Error: Cell angles have units that are not 'degrees' (%s)\n", angleUnits.c_str());
    return 1;
  }
  // Check units for lengths
  std::string lengthUnits = NC::GetAttrText(ncid_, cellLengthVID_, "units");
  if (setLengthFac(convert_h5_to_cpptraj_box_, lengthUnits, "Cell lengths"))
    return 1;
  // Get box lengths and angles to determine box type.
  start_[0] = 0;
  start_[1] = 0;
  count_[0] = 1; // 1 frame
  count_[1] = 3; // 3 coordinates (abg or xyz)
  double boxCrd[6]; /// XYZ ABG
  float* fptr = &ftmp_[0];
  if ( NC::CheckErr(nc_get_vara_float(ncid_, cellLengthVID_, start_, count_, fptr )) )
  {
    mprinterr("Error: Getting cell lengths.\n");
    return 1;
  }
  if ( NC::CheckErr(nc_get_vara_float(ncid_, cellAngleVID_, start_, count_, fptr+3)) )
  {
    mprinterr("Error: Getting cell angles.\n");
    return 1;
  }
  for (int i = 0; i < 6; i++)
    boxCrd[i] = (double)ftmp_[i];
  // Convert
  boxCrd[0] *= convert_h5_to_cpptraj_box_;
  boxCrd[1] *= convert_h5_to_cpptraj_box_;
  boxCrd[2] *= convert_h5_to_cpptraj_box_;
  mprintf("DEBUG:\tH5 Box: XYZ={%f %f %f} ABG={%f %f %f}\n",
           boxCrd[0], boxCrd[1], boxCrd[2], boxCrd[3], boxCrd[4], boxCrd[5]);
  if (ncbox.SetupFromXyzAbg( boxCrd )) {
    mprintf("Warning: H5 file unit cell variables appear to be empty; disabling box.\n");
    cellLengthVID_ = -1;
    cellAngleVID_ = -1;
    ncbox.SetNoBox();
    return -1;
  }
  // TODO check dim IDs

  return 0;
}

#endif /* HAS_HDF5 */


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

  // Set up coordinates
  int frameDID, atomDID, spatialDID, nframes;
  if (setupCoordVID(frameDID, atomDID, spatialDID, nframes)) {
    mprinterr("Error: Could not set up coordinates variable.\n");
    return TRAJIN_ERR;
  }
  mprintf("DEBUG: Unlimited dimid is %i\n", frameDID);
  mprintf("DEBUG: Atom dim is %i\n", atomDID);
  mprintf("DEBUG: Spatial dim is %i\n", spatialDID);

  // Check # atoms
  if ( natom_ != trajParm->Natom() ) {
    mprinterr("Error: Atom mismatch between topology (%i) and trajectory (%i).\n",
               trajParm->Natom(), natom_);
    return TRAJIN_ERR;
  }

  // Allocate temp space
  ftmp_.assign( natom_, 0 );

  // Check for box
  Box ncbox;
  int err = setupBoxVIDs(ncbox, frameDID, spatialDID);
  if (err == 1) {
    mprinterr("Error: Problem setting up box info.\n");
    return TRAJIN_ERR;
  }

  // Check for time
  err = nc_inq_varid(ncid_, "time", &timeVID_);
  if (err != NC_NOERR)
    timeVID_ = -1;

  // Get title
  SetTitle( NC::GetAttrText(ncid_, "TITLE") );

  // Setup coordinfo
  SetCoordInfo( CoordinateInfo( ncbox, false, false, (timeVID_ != -1) ) );

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
