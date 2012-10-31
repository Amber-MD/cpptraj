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
{
  //fprintf(stderr,"Amber Netcdf Restart Constructor\n");
  // Netcdf restarts always have 1 frame so always seekable
  seekable_=true;
}

// DESTRUCTOR
Traj_AmberRestartNC::~Traj_AmberRestartNC() {
  //fprintf(stderr,"Amber Netcdf Restart Destructor\n");
  // NOTE: Need to close file?
}

bool Traj_AmberRestartNC::ID_TrajFormat(CpptrajFile& fileIn) {
  if ( GetNetcdfConventions( FullFileStr() ) == NC_AMBERRESTART ) return true;
  return false;
}

// Traj_AmberRestartNC::closeTraj()
void Traj_AmberRestartNC::closeTraj() {
  NC_close();
}

// Traj_AmberRestartNC::openTraj()
/** Open up Netcdf restart file and set ncid. Variable and Dimension IDs are set 
  * up by SetupRead / SetupWrite and will not change for a given file between
  * open and close calls. 
  */
int Traj_AmberRestartNC::openTraj() {
  //mprintf("DEBUG: AmberRestartNC::open() called for %s, ncid=%i\n",BaseFileStr(),ncid);
  // If already open, return
  if (Ncid()!=-1) return 0;

  switch (access_) {

    case READ :
      if ( NC_openRead( FullFileStr() ) != 0 ) {
        mprinterr("Error: Opening Netcdf restart file %s for reading.\n",BaseFileStr());
        return 1;
      }
      break;
    
    case APPEND: 
      mprintf("Error: %s - Append is not supported by netcdf restart.\n",BaseFileStr());
      return 1;
      break;
    case WRITE:
      // All opening and closing of writes are done per-frame
      return 0;
  }

  if (debug_>0) mprintf("Successfully opened restart %s, ncid=%i\n",BaseFileStr(),Ncid());
  if (debug_>1) NetcdfDebug();

  return 0;
}

// Traj_AmberRestartNC::setupTrajin()
/** Set up netcdf restart file for reading, get all variable and dimension IDs. 
  * Also check number of atoms against associated parmtop.
  */
// NOTE: Replace attrText allocs with static buffer? 
int Traj_AmberRestartNC::setupTrajin(std::string const& fname, Topology* trajParm,
                    TrajInfo& tinfo)
{
  if (openTraj()) return -1;

  // Sanity check - Make sure this is a Netcdf restart
  if ( GetNetcdfConventions() != NC_AMBERRESTART ) {
    mprinterr("Error: Netcdf restart file %s conventions do not include \"AMBERRESTART\"\n",
            BaseFileStr());
    return -1;
  }

  // Get global attributes
  std::string attrText = GetAttrText("ConventionVersion");
  if (attrText!="1.0")
    mprintf("Warning: Netcdf restart file %s has ConventionVersion that is not 1.0 (%s)\n",
            BaseFileStr(), attrText.c_str());

  // Get title
  if (title_.empty())
    title_ = GetAttrText("title");

  // Setup Coordinates
  if ( SetupCoordinates()!=0 ) return -1;
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF restart file %s (%i) does not\n",
              BaseFileStr(),Ncatom());
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->Natom());
    return -1;
  }
  
  // Setup Velocity - allowed to fail gracefully
  if (SetupVelocity() == 0)
    hasVelocity_ = true;
  
  // Setup Time
  if ( SetupTime()!=0 ) return -1;
  
  // Box info
  int boxerr = SetupBox(boxAngle_, boxLength_);
  if (boxerr == 1)
    return 1;
  else if (boxerr == 0) {
    hasBox_ = true;
    // Get box lengths and angle so they can be checked
    count_[0] = 3;
    count_[1] = 0;
    if ( checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, boxLength_)) ) {
      mprinterr("Error: Getting cell lengths.\n"); 
      return 1;
    }
    if ( checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, boxAngle_)) ) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
 
  } else
    hasBox_ = false;

  // Replica Temperatures - allowed to fail gracefully
  if ( SetupTemperature() == 0 )
    hasTemperature_ = true;

  if ( SetupMultiD() == -1 ) return 1;

  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;

  closeTraj();
  return 1;
}

// Traj_AmberRestartNC::SetNoVelocity()
void Traj_AmberRestartNC::SetNoVelocity() {
  hasVelocity_=false;
}

// Traj_AmberRestartNC::processWriteArgs()
int processWriteArgs(ArgList& argIn) {
  // For write, assume we want velocities unless specified
  hasVelocity_=true;
  if (argIn->hasKey("novelocity")) this->SetNoVelocity();
  time0_ = argIn->getKeyDouble("time0", OUTPUTFRAMESHIFT);
  if (argIn->hasKey("remdtraj"))   this->SetTemperature();
  dt_ = argIn->getKeyDouble("dt",1.0);
  return 0;
}

// Traj_AmberRestartNC::setupTrajout()
/** Setting up is done for each frame.  */
int Traj_AmberRestartNC::setupTrajout(std::string const& fname, Topology* trajParm,
                     int NframesToWrite, TrajInfo const& tinfo)
{
  SetNcatom( trajParm->Natom() );
  //ncatom3 = ncatom * 3;
  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (NframesToWrite==1) singleWrite_=true;
  return 0;
}

// Traj_AmberRestartNC::setupWriteForSet()
/** Create Netcdf restart file for the given frame and set it up. Only set
  * up velocity info if both hasVelocity and V is not NULL.
  */
int Traj_AmberRestartNC::setupWriteForSet(int set, double *Vin) {
  std::string fname;
  // Create filename for this set
  // If just writing 1 frame dont modify output filename
  if (singleWrite_)
    fname = FullFileName();
  else
    fname = NumberFilename(FullFileName(), set+OUTPUTFRAMESHIFT);

  // Set up title
  if (title_.empty()) 
    title_.assign("Cpptraj Generated Restart");

  if ( NC_create( fname.c_str(), NC_AMBERRESTART, Ncatom(), (hasVelocity_ && Vin!=NULL),
                  hasBox_, hasTemperature_, (time0_>=0), title_ ) )
    return 1;

  // NOTE: Do not close here. Since this is called for every frame a write immediately
  //       follows. Close in writeFrame.
  
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
    if (debug_>1) mprintf("DEBUG: %s: Replica Temperature %lf\n",BaseFileStr(),T);
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
  if (hasVelocity_ && V!=NULL) {
    if ( checkNCerr(nc_get_vara_double(ncid_, velocityVID_, start_, count_, V)) ) {
      mprinterr("Error: Getting velocities\n"); 
      return 1;
    }
  }

  // Read box info 
  if (hasBox_) {
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
  if ( this->setupWriteForSet(set,V) ) return 1;

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
  if (hasVelocity_ && V!=NULL) {
    mprintf("DEBUG: Writing V, VID=%i\n",velocityVID_);
    if (checkNCerr(nc_put_vara_double(ncid_,velocityVID_,start_,count_,V)) ) {
      mprinterr("Error: Netcdf restart writing velocity %i\n",set);
      return 1;
    }
  }

  // write box
  if (hasBox_) { // && cellLengthVID!=-1) {
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
  if (time0_>=0) {
    restartTime_ = time0_;
    restartTime_ += (double) set;
    restartTime_ *= dt_;
    if (checkNCerr(nc_put_var_double(ncid_,timeVID_,&restartTime_)) ) {
      mprinterr("Error: Writing restart time.\n");
      return 1;
    }
  }

  // write temperature
  if (hasTemperature_) {
    if (checkNCerr(nc_put_var_double(ncid_,TempVID_,&T)) ) {
      mprinterr("Error: Writing restart temperature.\n"); 
      return 1;
    }
  }
  
  nc_sync(ncid_); // Necessary? File about to close anyway... 

  // Close file for this set
  closeTraj();

  return 0;
}  

// Traj_AmberRestartNC::info()
void Traj_AmberRestartNC::info() {
  mprintf("is a NetCDF AMBER restart file");
  if (hasVelocity_) mprintf(", with velocities");
  if (hasTemperature_) mprintf(", with replica temperature");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);

  /*if (debug_ > 2) {
      if (!title_.empt() )
        printfone("    title:        \"%s\"\n", title_.c_str());
      if (application != NULL)  
        printfone("    application:  \"%s\"\n", p->application);
      if (program != NULL) 
        printfone("    program:      \"%s\"\n", p->program);
      if (version != NULL) 
        printfone("    version:      \"%s\"\n", p->version);
  }*/
}
#endif
