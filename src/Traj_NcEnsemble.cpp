#ifdef BINTRAJ
#include "netcdf.h"
#include "Traj_NcEnsemble.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_NcEnsemble::Traj_NcEnsemble() : 
  Coord_(0),
  ensembleSize_(0),
  readAccess_(false),
  useVelAsCoords_(false)
{}

// DESTRUCTOR
Traj_NcEnsemble::~Traj_NcEnsemble() {
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
}

void Traj_NcEnsemble::WriteHelp() {
  mprintf("\t[remdtraj] [velocity]\n");
}

void Traj_NcEnsemble::ReadHelp() {
  mprintf("\t[usevelascoords]\n");
}

// Traj_NcEnsemble::processWriteArgs()
int Traj_NcEnsemble::processWriteArgs(ArgList& argIn) {
  SetTemperature(argIn.hasKey("remdtraj"));
  SetVelocity(argIn.hasKey("velocity"));
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
  if ( NC_openRead( filename_.Full() ) != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for reading.\n", filename_.base());
    return 1;
  }
  return 0;
}

// Traj_NcEnsemble::closeTraj()
void Traj_NcEnsemble::closeTraj() { NC_close(); }

// Traj_NcEnsemble::setupTrajin()
int Traj_NcEnsemble::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (filename_.SetFileNameWithExpansion( fname )) return TRAJIN_ERR;
  if (openTrajin()) return TRAJIN_ERR;
  readAccess_ = true;
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
  if ( SetupFrame()!=0 ) return TRAJIN_ERR;
  // Get ensemble info
  ensembleSize_ = SetupEnsembleDim();
  if (ensembleSize_ < 1) {
    mprinterr("Error: Could not get ensemble dimension info.\n");
    return TRAJIN_ERR;
  }
  // Setup Coordinates/Velocities
  if ( SetupCoordsVelo( useVelAsCoords_ )!=0 ) return TRAJIN_ERR;
  SetVelocity( HasVelocities() );
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF file %s (%i) does not\n"
              "Error:   match number in associated parmtop (%i)!\n",
              filename_.base(), Ncatom(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Setup Time
  if ( SetupTime()!=0 ) return TRAJIN_ERR;
  // Box info
  double boxcrd[6];
  if (SetupBox(boxcrd, NC_AMBERENSEMBLE) == 1) // 1 indicates an error
    return TRAJIN_ERR;
  SetBox( boxcrd );
  // Replica Temperatures - Allowed to fail silently
  if (SetupTemperature() == 0)
    SetTemperature( true );
  // Replica Dimensions
  ReplicaDimArray remdDim;
  if ( SetupMultiD(remdDim) == -1 ) return TRAJIN_ERR;
  SetReplicaDims( remdDim );
  // Allocate float array
  if (Coord_ != 0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  if (debug_>1) NetcdfDebug();
  closeTraj();
  return Ncframe();
}

// Traj_NcEnsemble::setupTrajout()
int Traj_NcEnsemble::setupTrajout(std::string const& fname, Topology* trajParm,
                                  int NframesToWrite, bool append)
{
  readAccess_ = false;
  if (!append) {
    ensembleSize_ = trajParm->EnsembleSize();
    filename_.SetFileName( fname );
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj Generated trajectory");
    // Create NetCDF file. TODO: Add option to set up replica indices.
    if ( NC_create( filename_.Full(), NC_AMBERENSEMBLE, trajParm->Natom(), HasV(),
                    false, HasBox(), HasT(), true, trajParm->ParmReplicaDimInfo(),
                    ensembleSize_, Title() ) )
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

// Traj_NcEnsemble::readArray()
int Traj_NcEnsemble::readArray(int set, FrameArray& f_ensemble) {
  start_[0] = set; // Frame
  start_[2] = 0;   // Atoms
  start_[3] = 0;   // XYZ
  count_[0] = 1;        // Frame
  count_[1] = 1;        // Ensemble
  count_[3] = 3;        // XYZ
  for (int member = 0; member != ensembleSize_; member++) {
    Frame& frm = f_ensemble[member];
    start_[1] = member;   // Ensemble
    count_[2] = Ncatom(); // Atoms
    // Read Coords
    if (checkNCerr(nc_get_vara_float(ncid_, coordVID_, start_, count_, Coord_)))
    {
      mprinterr("Error: Getting coordinates for frame %i\n", set+1);
      return 1;
    }
    FloatToDouble(frm.xAddress(), Coord_);
    //mprintf("Frm=%8i Rep=%8i ", set+1, member+1); // DEBUG
    //frm.printAtomCoord(0); // DEBUG
    // Read Velocities
    if (velocityVID_ != -1) {
      if (checkNCerr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Coord_)))
      {
        mprinterr("Error: Getting velocities for frame %i\n", set+1);
        return 1;
      }
      FloatToDouble(frm.vAddress(), Coord_);
    }
    // Read Box
    if (cellLengthVID_ != -1) {
      count_[2] = 3;
      if (checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, frm.bAddress())))
      {
        mprinterr("Error: Getting cell lengths for frame %i.\n", set+1);
        return 1;
      }
      if (checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, frm.bAddress()+3)))
      {
        mprinterr("Error: Getting cell angles for frame %i.\n", set+1);
        return 1;
      }
    }
    // Read Temperature
    if (TempVID_!=-1) {
      if (checkNCerr(nc_get_vara_double(ncid_, TempVID_, start_, count_, frm.tAddress())))
      {
        mprinterr("Error: Getting replica temperature for frame %i.\n", set+1);
        return 1;
      }
      //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
    }
    // Read indices
    if (indicesVID_!=-1) {
      count_[2] = remd_dimension_;
      if (checkNCerr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, frm.iAddress())))
      {
        mprinterr("Error: Getting replica indices for frame %i.\n", set+1);
        return 1;
      }
      //mprintf("DEBUG:\tReplica indices:");
      //for (int dim=0; dim < remd_dimension_; dim++) mprintf(" %i",remd_indices[dim]);
      //mprintf("\n");
    }
  }
  return 0;
}

// Traj_NcEnsemble::writeArray()
int Traj_NcEnsemble::writeArray(int set, FramePtrArray const& Farray) {
  start_[0] = ncframe_; // Frame
  start_[2] = 0;        // Atoms
  start_[3] = 0;        // XYZ
  count_[0] = 1; // Frame
  count_[1] = 1; // Ensemble
  count_[3] = 3; // XYZ
  for (int member = 0; member != ensembleSize_; member++) {
    Frame* frm = Farray[member];
    start_[1] = member;   // Ensemble
    count_[2] = Ncatom(); // Atoms
    // Write Coords
    //WriteIndices(); // DEBUG
    DoubleToFloat(Coord_, frm->xAddress());
    if (checkNCerr(nc_put_vara_float(ncid_, coordVID_, start_, count_, Coord_)))
    {
      mprinterr("Error: Netcdf Writing coords frame %i\n", set+1);
      return 1;
    }
    // Write velocity.
    if (HasV() && frm->HasVelocity()) { // TODO: Determine V beforehand
      DoubleToFloat(Coord_, frm->vAddress());
      if (checkNCerr(nc_put_vara_float(ncid_, velocityVID_, start_, count_, Coord_)) )
      {
        mprinterr("Error: Netcdf writing velocity frame %i\n", set+1);
        return 1;
      }
    }
    // Write box
    if (cellLengthVID_ != -1) {
      count_[2] = 3;
      if (checkNCerr(nc_put_vara_double(ncid_,cellLengthVID_,start_,count_,frm->bAddress())) )
      {
        mprinterr("Error: Writing cell lengths frame %i.\n", set+1);
        return 1;
      }
      if (checkNCerr(nc_put_vara_double(ncid_,cellAngleVID_,start_,count_,frm->bAddress()+3)))
      {
        mprinterr("Error: Writing cell angles frame %i.\n", set+1);
        return 1;
      }
    }
    // Write temperature
    if (TempVID_!=-1) {
      if (checkNCerr(nc_put_vara_double(ncid_,TempVID_,start_,count_,frm->tAddress())))
      {
        mprinterr("Error: Writing temperature frame %i.\n", set+1);
        return 1;
      }
    }
    // Write indices
    if (indicesVID_ != -1) {
      count_[2] = remd_dimension_;
      if (checkNCerr(nc_put_vara_int(ncid_,indicesVID_,start_,count_,frm->iAddress())))
      {
        mprinterr("Error: Writing indices frame %i.\n", set+1);
        return 1;
      }
    }
  }

  nc_sync(ncid_); // Necessary after every write??
  ++ncframe_;

  return 0;
}

// Traj_NcEnsemble::Info()
void Traj_NcEnsemble::Info() {
  mprintf("is a NetCDF Ensemble AMBER trajectory");
  if (readAccess_ && !HasCoords()) mprintf(" (no coordinates)");
  if (HasV()) mprintf(" containing velocities");
  if (HasT()) mprintf(" with replica temperatures");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);
}
#endif
