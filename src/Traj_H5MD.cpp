#include "Traj_H5MD.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Topology.h"
#ifdef BINTRAJ
# include <netcdf.h>
# include "NC_Routines.h"
# include "Units.h"
#endif
//#ifdef HAS_HDF5
//# include <H5Cpp.h>
//  using namespace H5;
//#endif

/// CONSTRUCTOR
Traj_H5MD::Traj_H5MD()
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
Traj_H5MD::~Traj_H5MD() {
  closeTraj();
//# ifdef HAS_HDF5
//  if (file_ != 0) delete file_;
//# endif
}


/** Identify trajectory format. File should be setup for READ */
bool Traj_H5MD::ID_TrajFormat(CpptrajFile& fileIn) {
# ifdef HAS_HDF5
  mprintf("BING BONG\n");
  int myNcid;
  if ( nc_open( fileIn.Filename().full(), NC_NOWRITE, &myNcid ) != NC_NOERR )
    return false;
  mainGroupNames_ = NC::GetGroupNames( myNcid, mainGroupIds_ );
  if (mainGroupNames_.empty()) return false;
  // Check for h5md and particle groups
  bool is_h5md = false;
  particle_gid_ = -1;
  for (unsigned int ii = 0; ii != mainGroupNames_.size(); ii++) {
    if (mainGroupNames_[ii] == "h5md") {
      is_h5md = true;
    } else if (mainGroupNames_[ii] == "particles") {
      particle_gid_ = mainGroupIds_[ii];
    }
  }
  nc_close( myNcid );
  if (is_h5md) {
    if (particle_gid_ == -1) {
      mprinterr("Error: H5MD file missing 'particle' group.\n");
    } else {
      return true;
    }
  }
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
void Traj_H5MD::Info() {
  mprintf("is a MDanalysis H5MD (HDF5) trajectory");
}

/** Close file. */
void Traj_H5MD::closeTraj() {
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
int Traj_H5MD::openTrajin() {

  return 0;
}

/** Read help */
void Traj_H5MD::ReadHelp() {

}

/** Process read arguments. */
int Traj_H5MD::processReadArgs(ArgList& argIn) {

  return 0;
}

# ifdef HAS_HDF5
/** Set up the coordinates variable ID, number of atoms, and number of frames. */
int Traj_H5MD::setupCoordVID(int position_gid, int& frameDID, int& atomDID,
                             int& spatialDID, int& nframes)
{
  // Get the 'coordinates' variable ID ('value')
  if (NC::CheckErr(nc_inq_varid(position_gid, "value", &coordVID_)))
    return 1;
  mprintf("DEBUG: Coordinates VID is %i\n", coordVID_);
  // Set conversion factor for coords
  std::string lengthUnits = NC::GetAttrText(position_gid, coordVID_, "unit");
  mprintf("DEBUG: length units are: %s\n", lengthUnits.c_str());
  if (Cpptraj::Units::SetConversionFactor( convert_h5_to_cpptraj_coord_, lengthUnits, "ang" )) {
    mprinterr("Error: Could not determine Coordinates conversion factor.\n");
    return 1;
  }

  // Dimensions
  frameDID = -1;
  atomDID = -1;
  spatialDID = -1;
  natom_ = 0;
  nframes = 0;
  // Need to get the unlimited dimension ID, which should be the frame dim.
  if (NC::CheckErr( nc_inq_unlimdim(position_gid, &frameDID) ) )
    return 1;
  if (frameDID < 0) {
    mprinterr("Error: No unlimited (frame) dimension present in H5MD file.\n");
    return 1;
  }
  // Get dimensions for coordinates
  int ndims = 0;
  if (NC::CheckErr( nc_inq_varndims(position_gid, coordVID_, &ndims) ) )
    return 1;
  if (ndims != 3) {
    mprinterr("Error: Expected 3 dims for 'coordinates', got %i\n", ndims);
    return 1;
  }
  int coord_dims[3];
  if (NC::CheckErr( nc_inq_vardimid(position_gid, coordVID_, coord_dims) ) )
    return 1;
  mprintf("DEBUG: Coord dims: %i %i %i\n", coord_dims[0], coord_dims[1], coord_dims[2]);
  // Check the dimensions. One should be frames (unlimited), one should be
  // atoms, and the last should be spatial (XYZ)

  bool dim_used[3];
  for (int i = 0; i < 3; i++) dim_used[i] = false;
  bool has_unlimited = false;
  for (int nd = 0; nd < ndims; nd++) {
    size_t dimsize = 0;
    if (NC::CheckErr( nc_inq_dimlen(position_gid, coord_dims[nd], &dimsize)))
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
  // TODO is this always the case for H5MD?
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
int Traj_H5MD::setupBoxVIDs(int box_gid, Box& ncbox, int frameDID, int spatialDID) {
  ncbox.SetNoBox();
  // Search for 'edges' group id
  int edges_gid = -1;
  Iarray box_ids;
  Sarray box_gnames = NC::GetGroupNames( box_gid, box_ids );
  for (unsigned int ii = 0; ii < box_gnames.size(); ii++) {
    if (box_gnames[ii] == "edges") {
      edges_gid = box_ids[ii];
    }
  }
  mprintf("DEBUG: edges group id %i\n", edges_gid);
  // Get the cell_lengths 'value' variable ID
  int err = nc_inq_varid(edges_gid, "value", &cellLengthVID_);
  if (err != NC_NOERR) return -1;
/*  // Get the 'cell_angles' variable ID 
  if (NC::CheckErr(nc_inq_varid(edges_gid, "cell_angles", &cellAngleVID_)))
    return 1;
  mprintf("DEBUG: Cell length vid= %i, cell angle vid= %i\n",
          cellLengthVID_, cellAngleVID_);
  // Ensure angles are in degrees
  std::string angleUnits = NC::GetAttrText(edges_gid, cellAngleVID_, "units");
  if (angleUnits != "degrees") {
    mprinterr("Error: Cell angles have units that are not 'degrees' (%s)\n", angleUnits.c_str());
    return 1;
  }*/
  // Check units for lengths
  std::string lengthUnits = NC::GetAttrText(edges_gid, cellLengthVID_, "unit");
  if (Cpptraj::Units::SetConversionFactor(convert_h5_to_cpptraj_box_, lengthUnits, "ang")) {
    mprinterr("Error: Could not determine Cell Lengths conversion factor.\n");
    return 1;
  }
  // Check what kind of info is stored
  int ndims = 0;
  if (NC::CheckErr( nc_inq_varndims(edges_gid, cellLengthVID_, &ndims) ) )
    return 1;
  if (ndims != 3) {
    mprinterr("Error: Expected 3 dims for 'edges', got %i\n", ndims);
    return 1;
  }

  // Get box lengths and angles to determine box type.
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1; // 1 frame
  count_[1] = 3; // 3 unit cell vectors
  count_[2] = 3; // 3 coordinates (XYZ)
  double ucell[9]; /// Unit cell vectors in rows
  float* fptr = &ftmp_[0];
  if ( NC::CheckErr(nc_get_vara_float(edges_gid, cellLengthVID_, start_, count_, fptr )) )
  {
    mprinterr("Error: Getting unit cell vectors.\n");
    return 1;
  }
  /*if ( NC::CheckErr(nc_get_vara_float(edges_gid, cellAngleVID_, start_, count_, fptr+3)) )
  {
    mprinterr("Error: Getting cell angles.\n");
    return 1;
  }*/
  // Convert
  for (int i = 0; i < 9; i++)
    ucell[i] = (double)ftmp_[i] * convert_h5_to_cpptraj_box_;
  // Convert
/*  boxCrd[0] *= convert_h5_to_cpptraj_box_;
  boxCrd[1] *= convert_h5_to_cpptraj_box_;
  boxCrd[2] *= convert_h5_to_cpptraj_box_;*/
  mprintf("DEBUG:\tH5MD Box: X={%f %f %f}\n"
          "      \t          Y={%f %f %f}\n"
          "      \t          Z={%f %f %f}\n",
           ucell[0], ucell[1], ucell[2],
           ucell[3], ucell[4], ucell[5],
           ucell[6], ucell[7], ucell[8]);
  if (ncbox.SetupFromUcell( ucell )) {
    mprintf("Warning: H5MD file unit cell variables appear to be empty; disabling box.\n");
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
int Traj_H5MD::setupTrajin(FileName const& fname, Topology* trajParm)
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

  // Get trajectory group in the particles group
  int trajectory_gid = -1;
  Iarray particle_ids;
  Sarray particle_gnames = NC::GetGroupNames( particle_gid_, particle_ids );
  for (unsigned int ii = 0; ii < particle_gnames.size(); ii++) {
    if (particle_gnames[ii] == "trajectory")
      trajectory_gid = particle_ids[ii];
  }
  if (trajectory_gid == -1) {
    mprinterr("Error: 'trajectory' group not found.\n");
    return TRAJIN_ERR;
  }
  // Get box and position groups in the trajectory group
  int box_gid = -1;
  int position_gid = -1;
  Iarray trajectory_ids;
  Sarray trajectory_gnames = NC::GetGroupNames( trajectory_gid, trajectory_ids );
  for (unsigned int ii = 0; ii < trajectory_gnames.size(); ii++) {
    if (trajectory_gnames[ii] == "box")
      box_gid = trajectory_ids[ii];
    else if (trajectory_gnames[ii] == "position")
      position_gid = trajectory_ids[ii];
  }
  mprintf("DEBUG: Box gid = %i, position gid = %i\n", box_gid, position_gid);

  // Set up coordinates
  int frameDID, atomDID, spatialDID, nframes;
  if (setupCoordVID(position_gid, frameDID, atomDID, spatialDID, nframes)) {
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
  ftmp_.assign( natom_*3, 0 );

  // Check for box
  Box ncbox;
  int err = setupBoxVIDs(box_gid, ncbox, frameDID, spatialDID);
  if (err == 1) {
    mprinterr("Error: Problem setting up box info.\n");
    return TRAJIN_ERR;
  }
  ncbox.PrintDebug("H5MD box");

  // Check for time
  err = nc_inq_varid(ncid_, "time", &timeVID_);
  if (err != NC_NOERR)
    timeVID_ = -1;

  // Get title
  SetTitle( NC::GetAttrText(ncid_, "TITLE") );

  // Setup coordinfo
  SetCoordInfo( CoordinateInfo( ncbox, false, false, (timeVID_ != -1) ) );

  return nframes;
# else
  return TRAJIN_ERR;
# endif
}

/** Read specified trajectory frame. */
int Traj_H5MD::readFrame(int set, Frame& frameIn) {
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = natom_;
  count_[2] = 3;

/*  // Get temperature
  if (TempVID_!=-1) {
    if ( NC::CheckErr(nc_get_vara_double(ncid_, TempVID_, start_, count_, frameIn.tAddress())) ) {
      mprinterr("Error: Getting replica temperature for frame %i.\n", set+1); 
      return 1;
    }
    //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
  }*/

  // Get time
  if (timeVID_!=-1) {
    float time;
    if (NC::CheckErr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time))) {
      mprinterr("Error: Getting time for frame %i.\n", set + 1);
      return 1;
    }
    frameIn.SetTime( (double)time );
  }

  float* fptr = &ftmp_[0];
  // Read Coords
  if (coordVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, coordVID_, start_, count_, fptr)) ) {
      mprinterr("Error: Getting coordinates for frame %i\n", set+1);
      return 1;
    }
    for (unsigned int idx = 0; idx != ftmp_.size(); idx++) 
      frameIn.xAddress()[idx] = (double)ftmp_[idx] * convert_h5_to_cpptraj_coord_;
  }

  // Read Velocities
  /*if (velocityVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting velocities for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.vAddress(), Coord_);
  }*/

  // Read Forces
  /*if (frcVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, frcVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting forces for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.fAddress(), Coord_);
  }*/

  // Read box info 
  if (cellLengthVID_ != -1) {
    double xyzabg[6];
    count_[1] = 3;
    count_[2] = 0;
    if (NC::CheckErr(nc_get_vara_float(ncid_, cellLengthVID_, start_, count_, fptr)))
    {
      mprinterr("Error: Getting cell lengths for frame %i.\n", set+1);
      return 1;
    }
    if (NC::CheckErr(nc_get_vara_float(ncid_, cellAngleVID_, start_, count_, fptr+3)))
    {
      mprinterr("Error: Getting cell angles for frame %i.\n", set+1);
      return 1;
    }
    for (int i = 0; i < 6; i++)
      xyzabg[i] = (double)ftmp_[i];
    // Convert
    xyzabg[0] *= convert_h5_to_cpptraj_box_;
    xyzabg[1] *= convert_h5_to_cpptraj_box_;
    xyzabg[2] *= convert_h5_to_cpptraj_box_;
    frameIn.ModifyBox().AssignFromXyzAbg( xyzabg );
  }

  return 0;
}

/** Read velocities from specified frame. */
int Traj_H5MD::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_H5MD::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_H5MD::WriteHelp() {

}

/** Process write arguments. */
int Traj_H5MD::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_H5MD::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{
  mprinterr("Error: H5 write not yet supported.\n");
  return 1;
}

/** Write specified trajectory frame. */
int Traj_H5MD::writeFrame(int set, Frame const& frameOut) {

  mprinterr("Error: H5 write not yet supported.\n");
  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_H5MD::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_H5MD::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_H5MD::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_H5MD::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_H5MD::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_H5MD::parallelCloseTraj() {

}
#endif
