#ifdef BINTRAJ
// This file contains a collection of routines designed for reading
// netcdf trajectory files used with amber.
// Dan Roe 10-2008
// Original implementation of netcdf in Amber by John Mongan.
#include "Traj_AmberNetcdf.h"
#include <netcdf.h>
#include "CpptrajStdio.h"
#include "NC_Routines.h"
#ifdef MPI
# include "ParallelNetcdf.h"
#endif

// CONSTRUCTOR
Traj_AmberNetcdf::Traj_AmberNetcdf() :
  Coord_(0),
  useVelAsCoords_(false),
  useFrcAsCoords_(false),
  readAccess_(false),
  outputTemp_(false),
  write_mdcrd_(false),
  write_mdvel_(false),
  write_mdfrc_(false)
{}

// DESTRUCTOR
Traj_AmberNetcdf::~Traj_AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
  // NOTE: Need to close file?
}

bool Traj_AmberNetcdf::ID_TrajFormat(CpptrajFile& fileIn) {
  if ( GetNetcdfConventions( fileIn.Filename().full() ) == NC_AMBERTRAJ ) return true;
  return false;
} 

// Traj_AmberNetcdf::close()
/** Close netcdf file. Set ncid to -1 since it can change between open
  * and close calls.
  */
void Traj_AmberNetcdf::closeTraj() {
  NC_close();
}

// Traj_AmberNetcdf::openTrajin()
int Traj_AmberNetcdf::openTrajin() {
  // If already open, return
  if (Ncid()!=-1) return 0;
  if ( NC_openRead( filename_.Full() ) != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for reading.\n", filename_.base()); 
    return 1;
  }
  return 0;
}

void Traj_AmberNetcdf::ReadHelp() {
  mprintf("\tusevelascoords: Use velocities instead of coordinates if present.\n"
          "\tusefrcascoords: Use forces instead of coordinates if present.\n");
}

int Traj_AmberNetcdf::processReadArgs(ArgList& argIn) {
  useVelAsCoords_ = argIn.hasKey("usevelascoords");
  useFrcAsCoords_ = argIn.hasKey("usefrcascoords");
  return 0;
}

// Traj_AmberNetcdf::setupTrajin()
/* * Open the netcdf file, read all dimension and variable IDs, close.
  * Return the number of frames in the file. 
  */
int Traj_AmberNetcdf::setupTrajin(FileName const& fname, Topology* trajParm)
{
  filename_ = fname;
  if (openTrajin()) return TRAJIN_ERR;
  readAccess_ = true;
  // Setup for Amber NetCDF trajectory
  if ( NC_setupRead(NC_AMBERTRAJ, trajParm->Natom(), useVelAsCoords_, useFrcAsCoords_) )
    return TRAJIN_ERR;
  // Get title
  SetTitle( GetNcTitle() );
  // Set coordinate info
  SetCoordInfo( NC_coordInfo() ); 
  // Amber Netcdf coords are float. Allocate a float array for converting
  // float to/from double.
  if (Coord_ != 0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  if (debug_>1) NetcdfDebug();
  closeTraj();
  return Ncframe();
}

void Traj_AmberNetcdf::WriteHelp() {
  mprintf("\tremdtraj: Write temperature to trajectory (makes REMD trajectory).\n"
          "\tmdvel   : Write only velocities to trajectory.\n"
          "\tmdfrc   : Write only forces to trajectory.\n"
          "\tmdcrd   : Write coordinates to trajectory (only required with mdvel/mdfrc).\n");
}

// Traj_AmberNetcdf::processWriteArgs()
int Traj_AmberNetcdf::processWriteArgs(ArgList& argIn) {
  outputTemp_ = argIn.hasKey("remdtraj");
  write_mdcrd_ = argIn.hasKey("mdcrd");
  if (argIn.hasKey("velocity"))
    mprintf("Warning: The 'velocity' keyword is no longer necessary and has been deprecated.\n");
  if (argIn.hasKey("force"))
    mprintf("Warning: The 'force' keyword is no longer necessary and has been deprecated.\n");
  write_mdvel_ = argIn.hasKey("mdvel");
  write_mdfrc_ = argIn.hasKey("mdfrc");
  return 0;
}

// Traj_AmberNetcdf::setupTrajout()
/** Create Netcdf file specified by filename and set up dimension and
  * variable IDs. 
  */
int Traj_AmberNetcdf::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{
  readAccess_ = false;
  if (!append) {
    CoordinateInfo cInfo = cInfoIn;
    // Deal with output options
    // For backwards compatibility always write temperature if remdtraj is true.
    if (outputTemp_ && !cInfo.HasTemp()) cInfo.SetTemperature(true);
    // Determine what kind of trajectory to write based on options/coordinate info.
    if (write_mdcrd_ || write_mdvel_ || write_mdfrc_) {
      // If an option was specified, check that CoordinateInfo has requested info.
      if (write_mdcrd_ && !cInfoIn.HasCrd()) { // SANITY CHECK
        mprinterr("Error: 'mdcrd' specified but no coordinate info present.\n");
        return 1;
      }
      if (write_mdvel_ && !cInfoIn.HasVel()) {
        mprinterr("Error: 'mdvel' specified but no velocity info present.\n");
        return 1;
      }
      if (write_mdfrc_ && !cInfoIn.HasForce()) {
        mprinterr("Error: 'mdfrc' specified but no force info present.\n");
        return 1;
      }
      cInfo.SetCrd( write_mdcrd_ );
      cInfo.SetVelocity( write_mdvel_ );
      cInfo.SetForce( write_mdfrc_ );
    }
    SetCoordInfo( cInfo );
    filename_ = fname;
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj Generated trajectory");
    // Create NetCDF file.
    if (NC_create( filename_.Full(), NC_AMBERTRAJ, trajParm->Natom(), CoordInfo(), Title() ))
      return 1;
    if (debug_>1) NetcdfDebug();
    // Close Netcdf file. It will be reopened write.
    NC_close();
    // Allocate memory
    if (Coord_!=0) delete[] Coord_;
    Coord_ = new float[ Ncatom3() ];
  } else { // NOTE: File existence is checked for in Trajout
    // Call setupTrajin to set input parameters. This will also allocate
    // memory for coords.
    if (setupTrajin(fname, trajParm) == TRAJIN_ERR) return 1;
    // Check output options.
    if (write_mdcrd_ || write_mdvel_ || write_mdfrc_)
      mprintf("Warning: 'mdcrd', 'mdvel', and 'mdfrc' are ignored for appending.\n");
    // Check that CoordinateInfo matches file we are appending to.
    if ((outputTemp_ || cInfoIn.HasTemp()) && !CoordInfo().HasTemp())
      mprintf("Warning: Cannot append temperature data to NetCDF file '%s'; no temperature dimension.\n",
              filename_.base());
    if (cInfoIn.HasVel() && !CoordInfo().HasVel())
      mprintf("Warning: Cannot append velocity data to NetCDF file '%s'; no velocity dimension.\n",
              filename_.base());
    if (cInfoIn.HasForce() && !CoordInfo().HasForce())
      mprintf("Warning: Cannot append force data to NetCDF file '%s'; no force dimension.\n",
              filename_.base());
    if (debug_ > 0)
      mprintf("\tNetCDF: Appending %s starting at frame %i\n", filename_.base(), Ncframe()); 
  }
  // Open file
  if ( NC_openWrite( filename_.Full() ) != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for Write.\n", filename_.base());
    return 1;
  }
  return 0;
}

// Traj_AmberNetcdf::readFrame()
/** Get the specified frame from amber netcdf file
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  */
int Traj_AmberNetcdf::readFrame(int set, Frame& frameIn) {
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;

  // Get temperature
  if (TempVID_!=-1) {
    if ( NC::CheckErr(nc_get_vara_double(ncid_, TempVID_, start_, count_, frameIn.tAddress())) ) {
      mprinterr("Error: Getting replica temperature for frame %i.\n", set+1); 
      return 1;
    }
    //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
  }

  // Get time
  if (timeVID_!=-1) {
    float time;
    if (NC::CheckErr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time))) {
      mprinterr("Error: Getting time for frame %i.\n", set + 1);
      return 1;
    }
    frameIn.SetTime( (double)time );
  }

  // Read Coords 
  if ( NC::CheckErr(nc_get_vara_float(ncid_, coordVID_, start_, count_, Coord_)) ) {
    mprinterr("Error: Getting coordinates for frame %i\n", set+1);
    return 1;
  }
  FloatToDouble(frameIn.xAddress(), Coord_);

  // Read Velocities
  if (velocityVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting velocities for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.vAddress(), Coord_);
  }

  // Read Forces
  if (frcVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, frcVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting forces for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.fAddress(), Coord_);
  }

  // Read indices. Input array must be allocated to be size remd_dimension.
  if (indicesVID_!=-1) {
    count_[1] = remd_dimension_;
    if ( NC::CheckErr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, frameIn.iAddress())) ) {
      mprinterr("Error: Getting replica indices for frame %i.\n", set+1);
      return 1;
    }
    //mprintf("DEBUG:\tReplica indices:");
    //for (int dim=0; dim < remd_dimension_; dim++) mprintf(" %i",remd_indices[dim]);
    //mprintf("\n");
  }

  // Read REMD values.
  ReadRemdValues(frameIn);

  // Read box info 
  if (cellLengthVID_ != -1) {
    count_[1] = 3;
    count_[2] = 0;
    if (NC::CheckErr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, frameIn.bAddress())))
    {
      mprinterr("Error: Getting cell lengths for frame %i.\n", set+1);
      return 1;
    }
    if (NC::CheckErr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, frameIn.bAddress()+3)))
    {
      mprinterr("Error: Getting cell angles for frame %i.\n", set+1);
      return 1;
    }
  }

  return 0;
}

// Traj_AmberNetcdf::readVelocity()
int Traj_AmberNetcdf::readVelocity(int set, Frame& frameIn) {
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  // Read Velocities
  if (velocityVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting velocities for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.vAddress(), Coord_);
  }
  return 0;
}

// Traj_AmberNetcdf::readForce()
int Traj_AmberNetcdf::readForce(int set, Frame& frameIn) {
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  // Read forces
  if (frcVID_ != -1) {
    if ( NC::CheckErr(nc_get_vara_float(ncid_, frcVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: Getting forces for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frameIn.fAddress(), Coord_);
  }
  return 0;
}

// Traj_AmberNetcdf::writeFrame() 
int Traj_AmberNetcdf::writeFrame(int set, Frame const& frameOut) {
  // Set indices for coords/velocities/forces
  start_[0] = ncframe_;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;

  // Write coords.
  if (coordVID_ != -1) {
    DoubleToFloat(Coord_, frameOut.xAddress());
    if (NC::CheckErr(nc_put_vara_float(ncid_,coordVID_,start_,count_,Coord_)) ) {
      mprinterr("Error: NetCDF writing coordinates frame %i\n", set+1);
      return 1;
    }
  }

  // Write velocity.
  if (velocityVID_ != -1) {
    DoubleToFloat(Coord_, frameOut.vAddress());
    if (NC::CheckErr(nc_put_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: NetCDF writing velocity frame %i\n", set+1);
      return 1;
    }
  }

  // Write forces.
  if (frcVID_ != -1) {
    DoubleToFloat(Coord_, frameOut.fAddress());
    if (NC::CheckErr(nc_put_vara_float(ncid_, frcVID_, start_, count_, Coord_)) ) {
      mprinterr("Error: NetCDF writing force frame %i\n", set+1);
      return 1;
    }
  }

  // Write Remd Values
  WriteRemdValues(frameOut);

  // Write box
  if (cellLengthVID_ != -1) {
    count_[1] = 3;
    count_[2] = 0;
    if (NC::CheckErr(nc_put_vara_double(ncid_, cellLengthVID_, start_, count_,
                                        frameOut.bAddress())) )
    {
      mprinterr("Error: Writing cell lengths frame %i.\n", set+1);
      return 1;
    }
    if (NC::CheckErr(nc_put_vara_double(ncid_, cellAngleVID_, start_, count_, 
                                        frameOut.bAddress()+3)) )
    {
      mprinterr("Error: Writing cell angles frame %i.\n", set+1);
      return 1;
    }
  }

  // Write temperature
  if (TempVID_ != -1) {
    if ( NC::CheckErr( nc_put_vara_double(ncid_,TempVID_,start_,count_,frameOut.tAddress())) ) {
      mprinterr("Error: Writing temperature frame %i.\n", set+1);
      return 1;
    }
  }

  // Write time
  if (timeVID_ != -1) {
    float tVal = (float)frameOut.Time();
    if ( NC::CheckErr( nc_put_vara_float(ncid_,timeVID_,start_,count_,&tVal)) ) {
      mprinterr("Error: Writing time frame %i.\n", set+1);
      return 1;
    }
  }
    
  // Write indices
  if (indicesVID_ != -1) {
    count_[1] = remd_dimension_;
    if ( NC::CheckErr(nc_put_vara_int(ncid_,indicesVID_,start_,count_,frameOut.iAddress())) ) {
      mprinterr("Error: Writing indices frame %i.\n", set+1);
      return 1;
    }
  }

  nc_sync(ncid_); // Necessary after every write??

  ++ncframe_;

  return 0;
}  

// Traj_AmberNetcdf::Info()
void Traj_AmberNetcdf::Info() {
  mprintf("is a NetCDF AMBER trajectory");
  if (readAccess_) {
    if (!HasCoords()) mprintf(" (no coordinates)");
    if (useVelAsCoords_) mprintf(" (using velocities as coordinates)");
    if (useFrcAsCoords_) mprintf(" (using forces as coordinates)");
    if (HasVelocities() || HasForces() || HasTemperatures()) {
      mprintf(" with");
      if (HasVelocities()) mprintf(" velocities");
      if (HasForces()) mprintf(" forces");
      if (HasTemperatures()) mprintf(" temperatures");
    }
    if (remd_dimension_ > 0) mprintf(" %i replica dimensions", remd_dimension_);
  } else {
    if ( !(write_mdcrd_ && write_mdvel_ && write_mdfrc_) ) {
      if (write_mdcrd_ || write_mdvel_ || write_mdfrc_) {
        mprintf(" with");
        if (write_mdcrd_) mprintf(" coordinates");
        if (write_mdvel_) mprintf(" velocities");
        if (write_mdfrc_) mprintf(" forces");
      }
    }
  }
}
#ifdef MPI
#ifdef HAS_PNETCDF
// =============================================================================
int Traj_AmberNetcdf::parallelOpenTrajin(Parallel::Comm const& commIn) {
  if (Ncid() != -1) return 0;
  int err = ncmpi_open(commIn.MPIcomm(), filename_.full(), NC_NOWRITE, MPI_INFO_NULL, &ncid_);
  if (checkPNCerr(err)) {
    mprinterr("Error: Opening NetCDF file %s for reading in parallel.\n", filename_.full());
    return 1;
  }
  err = ncmpi_begin_indep_data( ncid_ ); // Independent data mode
  return 0;
}

int Traj_AmberNetcdf::parallelOpenTrajout(Parallel::Comm const& commIn) {
  if (Ncid() != -1) return 0;
  int err = ncmpi_open(commIn.MPIcomm(), filename_.full(), NC_WRITE, MPI_INFO_NULL, &ncid_);
  if (checkPNCerr(err)) {
    mprinterr("Error: Opening NetCDF file '%s' for writing in parallel.\n", filename_.full());
    return 1;
  }
  err = ncmpi_begin_indep_data( ncid_ ); // Independent data mode
  return 0;
}

/** First master performs all necessary setup, then sends info to all children.
  */
int Traj_AmberNetcdf::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{
  int err = 0;
  if (commIn.Master()) {
    err = setupTrajout(fname, trajParm, cInfoIn, NframesToWrite, append);
    // NOTE: setupTrajout leaves file open. Should this change?
    NC_close();
  }
  commIn.MasterBcast(&err, 1, MPI_INT);
  if (err != 0) return 1;
  // Synchronize netcdf info on non-master threads.
  Sync(commIn);
  if (!commIn.Master()) {
    // Non masters need filename and allocate Coord
    filename_ = fname;
    if (Coord_ != 0) delete[] Coord_;
    Coord_ = new float[ Ncatom3() ];
  }
  return 0;
}

int Traj_AmberNetcdf::parallelReadFrame(int set, Frame& frameIn) {
  MPI_Offset pstart_[3];
  MPI_Offset pcount_[3];
  pstart_[0] = set;
  pstart_[1] = 0;
  pstart_[2] = 0;
  pcount_[0] = 1;
  pcount_[1] = Ncatom();
  pcount_[2] = 3;

  //int err = ncmpi_get_vara_float_all(ncid_, coordVID_, pstart_, pcount_, Coord_);
  int err = ncmpi_get_vara_float(ncid_, coordVID_, pstart_, pcount_, Coord_);
  if (checkPNCerr(err)) return Parallel::Abort(err);
  FloatToDouble(frameIn.xAddress(), Coord_);
  if (velocityVID_ != -1) {
    //err = ncmpi_get_vara_float_all(ncid_, velocityVID_, pstart_, pcount_, Coord_);
    err = ncmpi_get_vara_float(ncid_, velocityVID_, pstart_, pcount_, Coord_);
    if (checkPNCerr(err)) return Parallel::Abort(err);
    FloatToDouble(frameIn.vAddress(), Coord_);
  }
  if (frcVID_ != -1) {
    err = ncmpi_get_vara_float(ncid_, frcVID_, pstart_, pcount_, Coord_);
    if (checkPNCerr(err)) return Parallel::Abort(err);
    FloatToDouble(frameIn.fAddress(), Coord_);
  } 

  pcount_[2] = 0;
  if (cellLengthVID_ != -1) {
    pcount_[1] = 3;
    //err = ncmpi_get_vara_double_all(ncid_, cellLengthVID_, pstart_, pcount_, frameIn.bAddress());
    err = ncmpi_get_vara_double(ncid_, cellLengthVID_, pstart_, pcount_, frameIn.bAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
    //err = ncmpi_get_vara_double_all(ncid_, cellAngleVID_, pstart_, pcount_, frameIn.bAddress()+3);
    err = ncmpi_get_vara_double(ncid_, cellAngleVID_, pstart_, pcount_, frameIn.bAddress()+3);
  }
  if (TempVID_ != -1) {
    //err = ncmpi_get_vara_double_all(ncid_, TempVID_, pstart_, pcount_, frameIn.tAddress());
    err = ncmpi_get_vara_double(ncid_, TempVID_, pstart_, pcount_, frameIn.tAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  if (timeVID_ != -1) {
    float time;
    err = ncmpi_get_vara_float(ncid_, timeVID_, pstart_, pcount_, &time);
    if (checkPNCerr(err)) return Parallel::Abort(err);
    frameIn.SetTime( (double)time );
  }
  if (indicesVID_ != -1) {
    pcount_[1] = remd_dimension_;
    //err = ncmpi_get_vara_int_all(ncid_, indicesVID_, pstart_, pcount_, frameIn.iAddress());
    err = ncmpi_get_vara_int(ncid_, indicesVID_, pstart_, pcount_, frameIn.iAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  return 0;
}

int Traj_AmberNetcdf::parallelWriteFrame(int set, Frame const& frameOut) {
  MPI_Offset pstart_[3];
  MPI_Offset pcount_[3];
  pstart_[0] = set;
  pstart_[1] = 0;
  pstart_[2] = 0;
  pcount_[0] = 1;
  pcount_[1] = Ncatom();
  pcount_[2] = 3;
  // TODO check error better
  DoubleToFloat(Coord_, frameOut.xAddress());
  //int err = ncmpi_put_vara_float_all(ncid_, coordVID_, pstart_, pcount_, Coord_);
  int err = ncmpi_put_vara_float(ncid_, coordVID_, pstart_, pcount_, Coord_);
  if (checkPNCerr(err)) return Parallel::Abort(err);
  if (velocityVID_ != -1) {
    DoubleToFloat(Coord_, frameOut.vAddress());
    //err = ncmpi_put_vara_float_all(ncid_, velocityVID_, pstart_, pcount_, Coord_);
    err = ncmpi_put_vara_float(ncid_, velocityVID_, pstart_, pcount_, Coord_);
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  if (frcVID_ != -1) {
    DoubleToFloat(Coord_, frameOut.fAddress());
    err = ncmpi_put_vara_float(ncid_, frcVID_, pstart_, pcount_, Coord_);
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }

  pcount_[2] = 0;
  if (cellLengthVID_ != -1) {
    pcount_[1] = 3;
    //err = ncmpi_put_vara_double_all(ncid_, cellLengthVID_, pstart_, pcount_, frameOut.bAddress());
    err = ncmpi_put_vara_double(ncid_, cellLengthVID_, pstart_, pcount_, frameOut.bAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
    //err = ncmpi_put_vara_double_all(ncid_, cellAngleVID_, pstart_, pcount_, frameOut.bAddress()+3);
    err = ncmpi_put_vara_double(ncid_, cellAngleVID_, pstart_, pcount_, frameOut.bAddress()+3);
  }
  if (TempVID_ != -1) {
    //err = ncmpi_put_vara_double_all(ncid_, TempVID_, pstart_, pcount_, frameOut.tAddress());
    err = ncmpi_put_vara_double(ncid_, TempVID_, pstart_, pcount_, frameOut.tAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  if (timeVID_ != -1) {
    float tVal = (float)frameOut.Time();
    err = ncmpi_put_vara_float(ncid_, timeVID_, pstart_, pcount_, &tVal);
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  if (indicesVID_ != -1) {
    pcount_[1] = remd_dimension_;
    //err = ncmpi_put_vara_int_all(ncid_, indicesVID_, pstart_, pcount_, frameOut.iAddress());
    err = ncmpi_put_vara_int(ncid_, indicesVID_, pstart_, pcount_, frameOut.iAddress());
    if (checkPNCerr(err)) return Parallel::Abort(err);
  }
  return 0;
}

void Traj_AmberNetcdf::parallelCloseTraj() {
  if (ncid_ == -1) return;
  ncmpi_close( ncid_ );
  ncid_ = -1;
}
#else /* HAS_PNETCDF */
int Traj_AmberNetcdf::parallelOpenTrajin(Parallel::Comm const& commIn) { return 1; }
int Traj_AmberNetcdf::parallelOpenTrajout(Parallel::Comm const& commIn) { return 1; } 
int Traj_AmberNetcdf::parallelReadFrame(int set, Frame& frameIn) { return 1; }
int Traj_AmberNetcdf::parallelWriteFrame(int set, Frame const& frameOut) { return 1; }
void Traj_AmberNetcdf::parallelCloseTraj() { return; }
int Traj_AmberNetcdf::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{ 
  mprinterr("Error: NetCDF single trajectory output in parallel requires Pnetcdf.\n");
  return 1;
}
#endif /* HAS_PNETCDF */
#endif /* MPI */
#endif /* BINTRAJ */
