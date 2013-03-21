#ifdef BINTRAJ
// This file contains a collection of routines designed for reading
// (and writing?) netcdf restart files used with amber.
// Dan Roe 2011-01-07
#include "Traj_AmberRestartNC.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename
#include "netcdf.h"

// CONSTRUCTOR
Traj_AmberRestartNC::Traj_AmberRestartNC() :
  restartTime_(0),
  singleWrite_(false),
  time0_(OUTPUTFRAMESHIFT),
  dt_(1.0)
{ }

// DESTRUCTOR
Traj_AmberRestartNC::~Traj_AmberRestartNC() {
  //fprintf(stderr,"Amber Netcdf Restart Destructor\n");
  // NOTE: Need to close file?
}

bool Traj_AmberRestartNC::ID_TrajFormat(CpptrajFile& fileIn) {
  if ( GetNetcdfConventions( fileIn.Filename().full() ) == NC_AMBERRESTART ) return true;
  return false;
}

// Traj_AmberRestartNC::closeTraj()
void Traj_AmberRestartNC::closeTraj() {
  NC_close();
}

// Traj_AmberRestartNC::openTrajin()
int Traj_AmberRestartNC::openTrajin() {
  // If already open, return
  if (Ncid()!=-1) return 0;
  if ( NC_openRead( filename_.Full() ) != 0 ) {
    mprinterr("Error: Opening Netcdf restart file %s for reading.\n", filename_.base());
    return 1;
  }
  if (debug_>1) NetcdfDebug();
  return 0;
}

// Traj_AmberRestartNC::setupTrajin()
/** Set up netcdf restart file for reading, get all variable and dimension IDs. 
  * Also check number of atoms against associated parmtop.
  */
int Traj_AmberRestartNC::setupTrajin(std::string const& fname, Topology* trajParm)
{
  filename_.SetFileNameWithExpansion( fname );
  if (openTrajin()) return TRAJIN_ERR;
  // Sanity check - Make sure this is a Netcdf restart
  if ( GetNetcdfConventions() != NC_AMBERRESTART ) {
    mprinterr("Error: Netcdf restart file %s conventions do not include \"AMBERRESTART\"\n",
              filename_.base());
    return TRAJIN_ERR;
  }
  // Get global attributes
  std::string attrText = GetAttrText("ConventionVersion");
  if (attrText!="1.0")
    mprintf("Warning: Netcdf restart file %s has ConventionVersion that is not 1.0 (%s)\n",
            filename_.base(), attrText.c_str());
  // Get title
  SetTitle( GetAttrText("title") );
  // Setup Coordinates/Velocities
  if ( SetupCoordsVelo()!=0 ) return TRAJIN_ERR;
  SetVelocity( HasVelocities() );
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF restart file %s (%i) does not\n",
              filename_.base(), Ncatom());
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Setup Time
  if ( SetupTime()!=0 ) return TRAJIN_ERR;
  // Box info
  double boxcrd[6];
  if (SetupBox(boxcrd, NC_AMBERRESTART) == 1) // 1 indicates an error
    return TRAJIN_ERR;
  SetBox( boxcrd );
  // Replica Temperatures - allowed to fail silently 
  if (SetupTemperature() == 0)
    SetTemperature( true );
  if ( SetupMultiD() == -1 ) return TRAJIN_ERR;
  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;
  closeTraj();
  // Only 1 frame for NC restarts
  return 1;
}

// Traj_AmberRestartNC::processWriteArgs()
int Traj_AmberRestartNC::processWriteArgs(ArgList& argIn) {
  // For write, assume we want velocities unless specified
  SetVelocity(!argIn.hasKey("novelocity"));
  SetTemperature(argIn.hasKey("remdtraj"));
  time0_ = argIn.getKeyDouble("time0", 1.0);
  dt_ = argIn.getKeyDouble("dt",1.0);
  return 0;
}

// Traj_AmberRestartNC::setupTrajout()
/** Setting up is done for each frame.  */
int Traj_AmberRestartNC::setupTrajout(std::string const& fname, Topology* trajParm,
                                      int NframesToWrite, bool append)
{
  if (append) {
    mprinterr("Error: 'append' not supported by NetCDF restart\n");
    return 1;
  }
  filename_.SetFileName( fname );
  SetNcatom( trajParm->Natom() );
  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (NframesToWrite == 1) singleWrite_ = true;
  // Set up title
  if (Title().empty())
    SetTitle("Cpptraj Generated Restart");
  return 0;
}

// Traj_AmberRestartNC::readFrame()
/** Get the specified frame from amber netcdf file
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  */
int Traj_AmberRestartNC::readFrame(int set,double *X, double *V,double *box, double *T) {
  // Get temperature
  if (TempVID_!=-1) {
    if ( checkNCerr(nc_get_var_double(ncid_, TempVID_, T)) ) {
      mprinterr("Error: Getting replica temperature.\n");
      return 1;
    }
    if (debug_>1) mprintf("DEBUG: %s: Replica Temperature %lf\n",filename_.base(), T);
  }

  // Read Coords 
  start_[0] = 0;
  start_[1] = 0;
  count_[0] = Ncatom();
  count_[1] = 3;
  if ( checkNCerr(nc_get_vara_double(ncid_, coordVID_, start_, count_, X)) ) {
    mprinterr("Error: Getting Coords\n");
    return 1;
  }

  // Read Velocity
  if (velocityVID_!=-1 && V!=0) {
    if ( checkNCerr(nc_get_vara_double(ncid_, velocityVID_, start_, count_, V)) ) {
      mprinterr("Error: Getting velocities\n"); 
      return 1;
    }
  }

  // Read box info 
  if (cellLengthVID_ != -1) {
    count_[0] = 3;
    count_[1] = 0;
    if ( checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, box)) ) {
      mprinterr("Error: Getting cell lengths.\n"); 
      return 1;
    }
    if ( checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, box+3)) ) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
  }

  return 0;
}

int Traj_AmberRestartNC::readIndices(int set, int* remd_indices) {
  if (indicesVID_!=-1) {
    start_[0] = 0;
    count_[1] = remd_dimension_;
    if ( checkNCerr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, remd_indices)) ) {
      mprinterr("Error: Getting replica indices.\n");
      return 1;
    }
    //mprintf("DEBUG:\tReplica Rst indices:");
    //for (int dim=0; dim < remd_dimension_; dim++) mprintf(" %i",remd_indices[dim]);
    //mprintf("\n");
  }
  return 0;
}

// Traj_AmberRestartNC::writeFrame() 
int Traj_AmberRestartNC::writeFrame(int set, double *X, double *V,double *box, double T) {
  // Set up file for this set
  bool V_present = (HasV() && V != 0);
  std::string fname;
  // Create filename for this set
  // If just writing 1 frame dont modify output filename
  if (singleWrite_)
    fname = filename_.Full();
  else
    fname = NumberFilename(filename_.Full(), set+1);
  if ( NC_create( fname.c_str(), NC_AMBERRESTART, Ncatom(), V_present,
                  HasBox(), HasT(), (time0_ >= 0), Title() ) )
    return 1;
  // write coords
  start_[0] = 0;
  start_[1] = 0;
  count_[0] = Ncatom(); 
  count_[1] = 3;
  if (checkNCerr(nc_put_vara_double(ncid_,coordVID_,start_,count_,X)) ) {
    mprinterr("Error: Netcdf restart Writing coordinates %i\n",set);
    return 1;
  }
  // write velocity
  if (V_present) {
    mprintf("DEBUG: Writing V, VID=%i\n",velocityVID_);
    if (checkNCerr(nc_put_vara_double(ncid_,velocityVID_,start_,count_,V)) ) {
      mprinterr("Error: Netcdf restart writing velocity %i\n",set);
      return 1;
    }
  }
  // write box
  if (cellLengthVID_ != -1) {
    count_[0] = 3;
    count_[1] = 0;
    if (checkNCerr(nc_put_vara_double(ncid_,cellLengthVID_,start_,count_,box)) ) {
      mprinterr("Error: Writing cell lengths.\n");
      return 1;
    }
    if (checkNCerr(nc_put_vara_double(ncid_,cellAngleVID_,start_,count_, box+3)) ) {
      mprinterr("Error: Writing cell angles.\n");
      return 1;
    }
  }
  // write time
  if (timeVID_ != -1) {
    restartTime_ = (time0_ + (double)set) * dt_;
    if (checkNCerr(nc_put_var_double(ncid_,timeVID_,&restartTime_)) ) {
      mprinterr("Error: Writing restart time.\n");
      return 1;
    }
  }
  // write temperature
  if (TempVID_ != -1) {
    if (checkNCerr(nc_put_var_double(ncid_,TempVID_,&T)) ) {
      mprinterr("Error: Writing restart temperature.\n"); 
      return 1;
    }
  }
  //nc_sync(ncid_); // Necessary? File about to close anyway... 
  // Close file for this set
  closeTraj();
  return 0;
}  

// Traj_AmberRestartNC::Info()
void Traj_AmberRestartNC::Info() {
  mprintf("is a NetCDF AMBER restart file");
  if (HasV()) mprintf(", with velocities");
  if (HasT()) mprintf(", with replica temperature");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);

  /*if (debug_ > 2) {
      if (!title_.empt() )
        printfone("    title:        \"%s\"\n", title_.c_str());
      if (application != 0)  
        printfone("    application:  \"%s\"\n", p->application);
      if (program != 0) 
        printfone("    program:      \"%s\"\n", p->program);
      if (version != 0) 
        printfone("    version:      \"%s\"\n", p->version);
  }*/
}
#endif
