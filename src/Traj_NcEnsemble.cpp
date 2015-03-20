#ifdef BINTRAJ
#ifdef ENABLE_SINGLE_ENSEMBLE
//#include <cstdio> // DEBUG
#include "netcdf.h"
#include "Traj_NcEnsemble.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "MpiRoutines.h"
  // MPI cannot be defined for C++ when including mpi.h since it is reserved for the MPI namespace
# undef MPI
# define USE_MPI
# ifdef HAS_PNETCDF
#   include "pnetcdf.h"
# endif
#endif

// CONSTRUCTOR
Traj_NcEnsemble::Traj_NcEnsemble() : 
  Coord_(0),
  ensembleStart_(0),
  ensembleEnd_(0),
  readAccess_(false),
  useVelAsCoords_(false)
{}

// DESTRUCTOR
Traj_NcEnsemble::~Traj_NcEnsemble() {
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
}

void Traj_NcEnsemble::WriteHelp() {
  //mprintf("\t[remdtraj] [velocity]\n");
}

void Traj_NcEnsemble::ReadHelp() {
  mprintf("\t[usevelascoords]\n");
}

// Traj_NcEnsemble::processWriteArgs()
int Traj_NcEnsemble::processWriteArgs(ArgList& argIn) {
  //SetTemperature(argIn.hasKey("remdtraj"));
  //SetVelocity(argIn.hasKey("velocity"));
  return 0;
}

// Traj_NcEnsemble::processReadArgs()
int Traj_NcEnsemble::processReadArgs(ArgList& argIn) {
  useVelAsCoords_ = argIn.hasKey("usevelascoords");
  return 0;
}

// Traj_NcEnsemble::ID_TrajFormat()
bool Traj_NcEnsemble::ID_TrajFormat(CpptrajFile& fileIn) {
  return ( GetNetcdfConventions( fileIn.Filename().full() ) == NC_AMBERENSEMBLE );
}

// Traj_NcEnsemble::openTrajin()
int Traj_NcEnsemble::openTrajin() {
  // If already open, return
  if (Ncid()!=-1) return 0;
# if HAS_PNETCDF
  int err = ncmpi_open(MPI_COMM_WORLD, filename_.full(), NC_NOWRITE, MPI_INFO_NULL, &ncid_);
  //err += ncmpi_begin_indep_data( ncid_ ); // Disable independent data mode
# else
  int err = NC_openRead( filename_.Full() );
# endif
  if ( err != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for reading.\n", filename_.base());
    return 1;
  }
  return 0;
}

// Traj_NcEnsemble::closeTraj()
void Traj_NcEnsemble::closeTraj() {
# ifdef HAS_PNETCDF
  if (ncid_ == -1) return;
  //ncmpi_end_indep_data( ncid_ ); // Disable independent data mode.
  ncmpi_close( ncid_ ); 
  ncid_ = -1;
# else 
  NC_close();
# endif
}

#ifdef USE_MPI
static int NoPnetcdf() {
# ifndef HAS_PNETCDF
  mprinterr("Error: Compiled without pnetcdf. Netcdf Ensemble requires pnetcdf in parallel.\n");
  return 1;
# else
  return 0;
# endif
}
#endif

// Traj_NcEnsemble::setupTrajin()
int Traj_NcEnsemble::setupTrajin(std::string const& fname, Topology* trajParm)
{
# ifdef USE_MPI
  if (NoPnetcdf()) return TRAJIN_ERR;
# endif
  readAccess_ = true;
  if (filename_.SetFileNameWithExpansion( fname )) return TRAJIN_ERR;
  //if (openTrajin()) return TRAJIN_ERR;
  // Open single thread for now
  if (NC_openRead( filename_.Full() )) return TRAJIN_ERR;
  // Sanity check - Make sure this is a Netcdf ensemble trajectory
  if ( GetNetcdfConventions() != NC_AMBERENSEMBLE ) {
    mprinterr("Error: Netcdf file %s conventions do not include \"AMBERENSEMBLE\"\n",
              filename_.base());
    return TRAJIN_ERR;
  }
  // Get global attributes
  std::string attrText = GetAttrText("ConventionVersion");
  if ( attrText != "1.0")
    mprintf("Warning: Netcdf file %s has ConventionVersion that is not 1.0 (%s)\n",
            filename_.base(), attrText.c_str());
  // Get title
  SetTitle( GetAttrText("title") );
  // Get Frame info
  if ( SetupFrameDim()!=0 ) return TRAJIN_ERR;
  if ( Ncframe() < 1 ) {
    mprinterr("Error: Netcdf file is empty.\n");
    return TRAJIN_ERR;
  }
  // Get ensemble info
  int ensembleSize = SetupEnsembleDim();
  if (ensembleSize < 1) {
    mprinterr("Error: Could not get ensemble dimension info.\n");
    return TRAJIN_ERR;
  }
  // Setup Coordinates/Velocities
  if ( SetupCoordsVelo( useVelAsCoords_ )!=0 ) return TRAJIN_ERR;
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF file %s (%i) does not\n"
              "Error:   match number in associated parmtop (%i)!\n",
              filename_.base(), Ncatom(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Setup Time - FIXME: Allowed to fail silently
  SetupTime();
  // Box info
  double boxcrd[6];
  if (SetupBox(boxcrd, NC_AMBERENSEMBLE) == 1) // 1 indicates an error
    return TRAJIN_ERR;
  // Replica Temperatures - FIXME: Allowed to fail silently
  SetupTemperature();
  // Replica Dimensions
  ReplicaDimArray remdDim;
  if ( SetupMultiD(remdDim) == -1 ) return TRAJIN_ERR;
  // Set traj info: FIXME - no forces yet
  SetCoordInfo( CoordinateInfo(ensembleSize, remdDim, Box(boxcrd), HasVelocities(),
                               HasTemperatures(), HasTimes(), false) ); 
  if (debug_>1) NetcdfDebug();
  //closeTraj();
  // Close single thread for now
  NC_close();
  // Set up local ensemble parameters
# ifdef USE_MPI
  ensembleStart_ = worldrank;
  ensembleEnd_ = worldrank + 1;
# else
  ensembleStart_ = 0;
  ensembleEnd_ = ensembleSize;
# endif
  // DEBUG: Print info for all ranks
  WriteVIDs();
  // Allocate float array
  if (Coord_ != 0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  return Ncframe();
}

// Traj_NcEnsemble::setupTrajout()
int Traj_NcEnsemble::setupTrajout(std::string const& fname, Topology* trajParm,
                                  CoordinateInfo const& cInfoIn,
                                  int NframesToWrite, bool append)
{
  int err = 0;
# ifdef USE_MPI
  if (NoPnetcdf()) return 1;
# endif
  readAccess_ = false;
  if (!append) {
    CoordinateInfo cInfo = cInfoIn;
    // TODO: File output modifications
    SetCoordInfo( cInfo );
#   ifdef USE_MPI
    ensembleStart_ = worldrank;
    ensembleEnd_ = worldrank + 1;
#   else
    ensembleStart_ = 0;
    ensembleEnd_ = cInfo.EnsembleSize();;
#   endif
    filename_.SetFileName( fname );
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj Generated trajectory");
#   ifdef USE_MPI
    if (worldrank == 0) { // Only master creates file.
#   endif
      // Create NetCDF file.
      err = NC_create(filename_.Full(), NC_AMBERENSEMBLE, trajParm->Natom(), CoordInfo(), Title());
      if (debug_ > 1 && err == 0) NetcdfDebug();
      // Close Netcdf file. It will be reopened write.
      NC_close();
#   ifdef USE_MPI
    }
    parallel_bcastMaster(&err, 1, PARA_INT);
#   endif
    if (err != 0) return 1;
#   ifdef USE_MPI
    // Synchronize netcdf info on non-master threads
    Sync();
    // DEBUG: Print info for all ranks
    WriteVIDs();
#   endif
    // Allocate memory
    if (Coord_!=0) delete[] Coord_;
    Coord_ = new float[ Ncatom3() ];
  } else { // NOTE: File existence is checked for in Trajout
    // Call setupTrajin to set input parameters. This will also allocate
    // memory for coords.
    if (setupTrajin(fname, trajParm) == TRAJIN_ERR) return 1;
    if (debug_ > 0)
      mprintf("\tNetCDF: Appending %s starting at frame %i\n", filename_.base(), Ncframe());
  }
  // Open file
# ifdef HAS_PNETCDF
  err = ncmpi_open(MPI_COMM_WORLD, filename_.full(), NC_WRITE, MPI_INFO_NULL, &ncid_);
  // TODO: Graceful error handling
# else
  err = NC_openWrite( filename_.Full() );
# endif  
  if ( err != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for Write.\n", filename_.base());
    return 1;
  }

  return 0;
}

// Traj_NcEnsemble::readFrame()
int Traj_NcEnsemble::readFrame(int set, Frame& frameIn) {
  return 1;
}

// Traj_NcEnsemble::readVelocity()
int Traj_NcEnsemble::readVelocity(int set, Frame& frameIn) {
  return 1;
}

// Traj_NcEnsemble::writeFrame()
int Traj_NcEnsemble::writeFrame(int set, Frame const& frameOut) {
  return 1;
}

#ifdef HAS_PNETCDF
inline int checkPNCerr(int err) {
  if (err != NC_NOERR) {
    rprinterr("PnetCDF Error: %s\n", ncmpi_strerror(err));
    parallel_abort( err );
    return 1;
  }
  return 0;
}
#endif

// Traj_NcEnsemble::readArray()
int Traj_NcEnsemble::readArray(int set, FrameArray& f_ensemble) {
# ifdef HAS_PNETCDF
  MPI_Offset pstart_[4];
  MPI_Offset pcount_[4];
# define start_ pstart_
# define count_ pcount_
# endif
  start_[0] = set; // Frame
  start_[2] = 0;   // Atoms
  start_[3] = 0;   // XYZ
  count_[0] = 1;        // Frame
  count_[1] = 1;        // Ensemble
  count_[3] = 3;        // XYZ
  //rprintf("DEBUG: Reading frame %i\n", set+1);
  for (int member = ensembleStart_; member != ensembleEnd_; member++) {
#   ifdef USE_MPI
    Frame& frm = f_ensemble[0];
#   else
    Frame& frm = f_ensemble[member];
#   endif
    start_[1] = member;   // Ensemble
    count_[2] = Ncatom(); // Atoms
    // Read Coords
#   ifdef HAS_PNETCDF
    if (checkPNCerr(ncmpi_get_vara_float_all(ncid_, coordVID_, start_, count_, Coord_)))
#   else
    if (checkNCerr(nc_get_vara_float(ncid_, coordVID_, start_, count_, Coord_)))
#   endif
    {
      rprinterr("Error: Getting coordinates for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frm.xAddress(), Coord_);
    //mprintf("Frm=%8i Rep=%8i ", set+1, member+1); // DEBUG
    //frm.printAtomCoord(0); // DEBUG
    // Read Velocities
    if (velocityVID_ != -1) {
#     ifdef HAS_PNETCDF
      if (checkPNCerr(ncmpi_get_vara_float_all(ncid_, velocityVID_, start_, count_, Coord_)))
#     else
      if (checkNCerr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Coord_)))
#     endif
      {
        rprinterr("Error: Getting velocities for frame %i\n", set+1);
        return 1;
      }
      FloatToDouble(frm.vAddress(), Coord_);
    }
    // Read Box
    if (cellLengthVID_ != -1) {
      count_[2] = 3;
#     ifdef HAS_PNETCDF
      if (checkPNCerr(ncmpi_get_vara_double_all(ncid_, cellLengthVID_, start_, count_, frm.bAddress())))
#     else
      if (checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, frm.bAddress())))
#     endif
      {
        rprinterr("Error: Getting cell lengths for frame %i.\n", set+1);
        return 1;
      }
#     ifdef HAS_PNETCDF
      if (checkPNCerr(ncmpi_get_vara_double_all(ncid_, cellAngleVID_, start_, count_, frm.bAddress()+3)))
#     else
      if (checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, frm.bAddress()+3)))
#     endif
      {
        rprinterr("Error: Getting cell angles for frame %i.\n", set+1);
        return 1;
      }
    }
    // Read Temperature
    if (TempVID_!=-1) {
#     ifdef HAS_PNETCDF
      if (checkPNCerr(ncmpi_get_vara_double_all(ncid_, TempVID_, start_, count_, frm.tAddress())))
#     else
      if (checkNCerr(nc_get_vara_double(ncid_, TempVID_, start_, count_, frm.tAddress())))
#     endif
      {
        rprinterr("Error: Getting replica temperature for frame %i.\n", set+1);
        return 1;
      }
      //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
    }
    // Read indices
    if (indicesVID_!=-1) {
      count_[2] = remd_dimension_;
#     ifdef HAS_PNETCDF
      if (checkPNCerr(ncmpi_get_vara_int_all(ncid_, indicesVID_, start_, count_, frm.iAddress())))
#     else
      if (checkNCerr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, frm.iAddress())))
#     endif
      {
        rprinterr("Error: Getting replica indices for frame %i.\n", set+1);
        return 1;
      }
      // DEBUG
      //char buffer[128];
      //char* ptr = buffer;
      //ptr += sprintf(buffer,"DEBUG:\tReplica indices:");
      //for (int dim=0; dim < remd_dimension_; dim++) ptr += sprintf(ptr, " %i", frm.RemdIndices()[dim]);
      //sprintf(ptr,"\n");
      //rprintf("%s", buffer);
    }
  }
# ifdef HAS_PNETCDF
  // DEBUG
# undef start_
# undef count_
# endif
  return 0;
}

// Traj_NcEnsemble::writeArray()
int Traj_NcEnsemble::writeArray(int set, FramePtrArray const& Farray) {
# ifdef HAS_PNETCDF
  MPI_Offset pstart_[4];
  MPI_Offset pcount_[4];
# define start_ pstart_
# define count_ pcount_
# endif
  start_[0] = ncframe_; // Frame
  start_[2] = 0;        // Atoms
  start_[3] = 0;        // XYZ
  count_[0] = 1; // Frame
  count_[1] = 1; // Ensemble
  count_[3] = 3; // XYZ
  for (int member = ensembleStart_; member != ensembleEnd_; member++) {
    //rprintf("DEBUG: Writing set %i, member %i\n", set+1, member); 
#   ifdef USE_MPI
    Frame* frm = Farray[0];
#   else
    Frame* frm = Farray[member];
#   endif
    start_[1] = member;   // Ensemble
    count_[2] = Ncatom(); // Atoms
    // Write Coords
    //WriteIndices(); // DEBUG
    DoubleToFloat(Coord_, frm->xAddress());
#   ifdef HAS_PNETCDF
    if (ncmpi_put_vara_float_all(ncid_, coordVID_, start_, count_, Coord_))
#   else
    if (checkNCerr(nc_put_vara_float(ncid_, coordVID_, start_, count_, Coord_)))
#   endif
    {
      mprinterr("Error: Netcdf Writing coords frame %i\n", set+1);
      return 1;
    }
    // Write velocity.
    if (CoordInfo().HasVel() && frm->HasVelocity()) { // TODO: Determine V beforehand
      DoubleToFloat(Coord_, frm->vAddress());
#     ifdef HAS_PNETCDF
      if (ncmpi_put_vara_float_all(ncid_, velocityVID_, start_, count_, Coord_))
#     else
      if (checkNCerr(nc_put_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) )
#     endif
      {
        mprinterr("Error: Netcdf writing velocity frame %i\n", set+1);
        return 1;
      }
    }
    // Write box
    if (cellLengthVID_ != -1) {
      count_[2] = 3;
#     ifdef HAS_PNETCDF
      if (ncmpi_put_vara_double_all(ncid_,cellLengthVID_,start_,count_,frm->bAddress()))
#     else
      if (checkNCerr(nc_put_vara_double(ncid_,cellLengthVID_,start_,count_,frm->bAddress())) )
#     endif
      {
        mprinterr("Error: Writing cell lengths frame %i.\n", set+1);
        return 1;
      }
#     ifdef HAS_PNETCDF
      if (ncmpi_put_vara_double_all(ncid_,cellAngleVID_,start_,count_,frm->bAddress()+3))
#     else
      if (checkNCerr(nc_put_vara_double(ncid_,cellAngleVID_,start_,count_,frm->bAddress()+3)))
#     endif
      {
        mprinterr("Error: Writing cell angles frame %i.\n", set+1);
        return 1;
      }
    }
    // Write temperature
    if (TempVID_!=-1) {
#     ifdef HAS_PNETCDF
      if (ncmpi_put_vara_double_all(ncid_,TempVID_,start_,count_,frm->tAddress()))
#     else
      if (checkNCerr(nc_put_vara_double(ncid_,TempVID_,start_,count_,frm->tAddress())))
#     endif
      {
        mprinterr("Error: Writing temperature frame %i.\n", set+1);
        return 1;
      }
    }
    // Write indices
    if (indicesVID_ != -1) {
      count_[2] = remd_dimension_;
#     ifdef HAS_PNETCDF
      if (ncmpi_put_vara_int_all(ncid_,indicesVID_,start_,count_,frm->iAddress()))
#     else
      if (checkNCerr(nc_put_vara_int(ncid_,indicesVID_,start_,count_,frm->iAddress())))
#     endif
      {
        mprinterr("Error: Writing indices frame %i.\n", set+1);
        return 1;
      }
    }
  }
# ifdef HAS_PNETCDF
  //ncmpi_sync(ncid_);
# else
  nc_sync(ncid_); // Necessary after every write??
# endif
  ++ncframe_;
# ifdef HAS_PNETCDF
  // DEBUG
# undef start_
# undef count_
# endif
  return 0;
}

// Traj_NcEnsemble::Info()
void Traj_NcEnsemble::Info() {
  mprintf("is a NetCDF Ensemble AMBER trajectory");
  if (readAccess_ && !HasCoords()) mprintf(" (no coordinates)");
  //if (HasV()) mprintf(" containing velocities");
  //if (HasT()) mprintf(" with replica temperatures");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);
}
#endif
#endif
