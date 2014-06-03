#ifdef BINTRAJ
#include "netcdf.h"
#include "Traj_NcEnsemble.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_NcEnsemble::Traj_NcEnsemble() : 
  Coord_(0),
  ensembleSize_(0),
  readAccess_(false) {}

// DESTRUCTOR
Traj_NcEnsemble::~Traj_NcEnsemble() {
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
}

bool Traj_NcEnsemble::ID_TrajFormat(CpptrajFile& fileIn) {
  return ( GetNetcdfConventions( fileIn.Filename().full() ) == NC_AMBERENSEMBLE );
}

void Traj_NcEnsemble::closeTraj() { NC_close(); }

int Traj_NcEnsemble::openTrajin() {
  // If already open, return
  if (Ncid()!=-1) return 0;
  if ( NC_openRead( filename_.Full() ) != 0 ) {
    mprinterr("Error: Opening Netcdf file %s for reading.\n", filename_.base());
    return 1;
  }
  return 0;
}

int Traj_NcEnsemble::setupTrajin(std::string const& fname, Topology* trajParm)
{
  return TRAJIN_ERR;
}

int Traj_NcEnsemble::processWriteArgs(ArgList& argIn) {
  return 0;
}

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
    //if (debug_>1)
      NetcdfDebug();
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

int Traj_NcEnsemble::readFrame(int set, Frame& frameIn) {
  return 1;
}

int Traj_NcEnsemble::readVelocity(int set, Frame& frameIn) {
  return 1;
}

int Traj_NcEnsemble::writeFrame(int set, Frame const& frameOut) {
  return 1;
}

int Traj_NcEnsemble::writeArray(int set, FramePtrArray const& Farray) {
  start_[1] = ncframe_;
  start_[2] = 0;
  start_[3] = 0;
  count_[0] = 1;        // Ensemble
  count_[1] = 1;        // Frame
  count_[2] = Ncatom(); // Atoms
  count_[3] = 3;        // XYZ
  for (int member = 0; member != ensembleSize_; member++) {
    start_[0] = member;
    // Coords
    DoubleToFloat(Coord_, Farray[member]->xAddress());
    if (checkNCerr(nc_put_vara_float(ncid_, coordVID_, start_, count_, Coord_))) {
      mprinterr("Error: Netcdf Writing coords frame %i\n", set+1);
      return 1;
    }
    // Write velocity.
    // Write box
    // Write temperature
    // Write indices

  }

  nc_sync(ncid_); // Necessary after every write??
  ++ncframe_;

  return 0;
}

void Traj_NcEnsemble::Info() {
  mprintf("is a NetCDF Ensemble AMBER trajectory");
  if (readAccess_ && !HasCoords()) mprintf(" (no coordinates)");
  if (HasV()) mprintf(" containing velocities");
  if (HasT()) mprintf(" with replica temperatures");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);
}
#endif
