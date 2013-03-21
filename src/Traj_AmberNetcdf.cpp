#ifdef BINTRAJ
// This file contains a collection of routines designed for reading
// netcdf trajectory files used with amber.
// Dan Roe 10-2008
// Original implementation of netcdf in Amber by John Mongan.
#include "netcdf.h"
#include "Traj_AmberNetcdf.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_AmberNetcdf::Traj_AmberNetcdf() :
  Coord_(0),
  Veloc_(0)
{ }

// DESTRUCTOR
Traj_AmberNetcdf::~Traj_AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
  if (Veloc_!=0) delete[] Veloc_;
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

// Traj_AmberNetcdf::setupTrajin()
/* * Open the netcdf file, read all dimension and variable IDs, close.
  * Return the number of frames in the file. 
  */
int Traj_AmberNetcdf::setupTrajin(std::string const& fname, Topology* trajParm)
{
  filename_.SetFileNameWithExpansion( fname );
  if (openTrajin()) return TRAJIN_ERR;

  // Sanity check - Make sure this is a Netcdf trajectory
  if ( GetNetcdfConventions() != NC_AMBERTRAJ ) {
    mprinterr("Error: Netcdf file %s conventions do not include \"AMBER\"\n",filename_.base());
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
  // Setup Coordinates/Velocities
  if ( SetupCoordsVelo()!=0 ) return TRAJIN_ERR;
  SetVelocity( HasVelocities() );
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF file %s (%i) does not\n",
              filename_.base(),Ncatom());
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Setup Time
  if ( SetupTime()!=0 ) return TRAJIN_ERR;
  // Box info
  double boxcrd[6];
  if (SetupBox(boxcrd, NC_AMBERTRAJ) == 1) // 1 indicates an error
    return TRAJIN_ERR;
  SetBox( boxcrd );
  // Replica Temperatures - Allowed to fail silently
  if (SetupTemperature() == 0)
    SetTemperature( true );
  // Replica Dimensions
  if ( SetupMultiD() == -1 ) return TRAJIN_ERR;
  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;
  // Amber Netcdf coords are float. Allocate a float array for converting
  // float to/from double.
  if (Coord_ != 0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  if (Veloc_ != 0) delete[] Veloc_;
  if (velocityVID_ != -1) 
    Veloc_ = new float[ Ncatom3() ];
  else
    Veloc_ = 0;
  if (debug_>1) NetcdfDebug();
  closeTraj();
  // NetCDF files are always seekable
  SetSeekable( true );
  return Ncframe();
}

// Traj_AmberNetcdf::processWriteArgs()
int Traj_AmberNetcdf::processWriteArgs(ArgList& argIn) {
  SetTemperature(argIn.hasKey("remdtraj"));
  SetVelocity(argIn.hasKey("velo"));
  return 0;
}

// Traj_AmberNetcdf::setupTrajout()
/** Create Netcdf file specified by filename and set up dimension and
  * variable IDs. 
  */
int Traj_AmberNetcdf::setupTrajout(std::string const& fname, Topology* trajParm, 
                                   int NframesToWrite, bool append)
{
  if (!append) {
    filename_.SetFileName( fname );
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj Generated trajectory");
    // Create NetCDF file.
    if ( NC_create( filename_.Full(), NC_AMBERTRAJ, trajParm->Natom(), HasV(),
                    HasBox(), HasT(), true, Title() ) )
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

// Traj_AmberNetcdf::readFrame()
/** Get the specified frame from amber netcdf file
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  */
int Traj_AmberNetcdf::readFrame(int set,double *X, double *V,double *box, double *T) {
  // Get temperature
  if (TempVID_!=-1) {
    start_[0] = set;
    count_[0] = 1;
    if ( checkNCerr(nc_get_vara_double(ncid_, TempVID_, start_, count_, T)) ) {
      mprinterr("Error: Getting replica temperature.\n"); 
      return 1;
    }
    //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
  }

  // Read Coords 
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  if ( checkNCerr(nc_get_vara_float(ncid_, coordVID_, start_, count_, Coord_)) ) {
    mprinterr("Error: Getting frame %i\n", set);
    return 1;
  }
  FloatToDouble(X, Coord_);

  // Read Velocities
  if (velocityVID_ != -1) {
    if ( checkNCerr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Veloc_)) ) {
      mprinterr("Error: Getting velocities for frame %i\n", set);
      return 1;
    }
    FloatToDouble(V, Veloc_);
  }

  // Read box info 
  if (cellLengthVID_ != -1) {
    count_[1] = 3;
    count_[2] = 0;
    if ( checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, box)) ) {
      mprinterr("Getting cell lengths.\n");
      return 1;
    }
    if ( checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, box+3)) ) {
      mprinterr("Getting cell angles.\n");
      return 1;
    }
  }


  return 0;
}

int Traj_AmberNetcdf::readVelocity(int set, double* V) {
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  // Read Velocities
  if (velocityVID_ != -1) {
    if ( checkNCerr(nc_get_vara_float(ncid_, velocityVID_, start_, count_, Veloc_)) ) {
      mprinterr("Error: Getting velocities for frame %i\n", set);
      return 1;
    }
    FloatToDouble(V, Veloc_);
  }
  return 0;
}

/** Read REMD indices. Input array must be allocated to be size remd_dimension
  * by prior call to NreplicaDimensions().
  */
int Traj_AmberNetcdf::readIndices(int set, int* remd_indices) {
  if (indicesVID_!=-1) {
    start_[0] = set;
    start_[1] = 0;
    count_[0] = 1;
    count_[1] = remd_dimension_;
    if ( checkNCerr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, remd_indices)) ) {
      mprinterr("Error: Getting replica indices.\n");
      return 1;
    }
    //mprintf("DEBUG:\tReplica indices:");
    //for (int dim=0; dim < remd_dimension_; dim++) mprintf(" %i",remd_indices[dim]);
    //mprintf("\n");
  }
  return 0;
}

// Traj_AmberNetcdf::writeFrame() 
int Traj_AmberNetcdf::writeFrame(int set, double *X, double *V, double *box, double T) {

  DoubleToFloat(Coord_, X);

  // Write coords
  start_[0] = ncframe_;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  if (checkNCerr(nc_put_vara_float(ncid_,coordVID_,start_,count_,Coord_)) ) {
    mprinterr("Error: Netcdf Writing frame %i\n",set);
    return 1;
  }

  // Write box
  if (cellLengthVID_ != -1) {
    count_[1] = 3;
    count_[2] = 0;
    if (checkNCerr(nc_put_vara_double(ncid_,cellLengthVID_,start_,count_,box)) ) {
      mprinterr("Error: Writing cell lengths.\n");
      return 1;
    }
    if (checkNCerr(nc_put_vara_double(ncid_,cellAngleVID_,start_,count_, box+3)) ) {
      mprinterr("Error: Writing cell angles.\n");
      return 1;
    }
  }

  // Write temperature
  if (TempVID_!=-1) {
    if ( checkNCerr( nc_put_vara_double(ncid_,TempVID_,start_,count_,&T)) ) {
      mprinterr("Error: Writing temperature.\n");
      return 1;
    }
  }
  
  nc_sync(ncid_); // Necessary after every write??

  ++ncframe_;

  return 0;
}  

// Traj_AmberNetcdf::info()
void Traj_AmberNetcdf::Info() {
  mprintf("is a NetCDF AMBER trajectory");
  if (HasV()) mprintf(" containing velocities");
  if (HasT()) mprintf(" with replica temperatures");
  if (remd_dimension_ > 0) mprintf(", with %i dimensions", remd_dimension_);

  /*if (debug_ > 2) {
      if (title != 0)
        printfone("    title:        \"%s\"\n", p->title);
      if (application != 0)  
        printfone("    application:  \"%s\"\n", p->application);
      if (program != 0) 
        printfone("    program:      \"%s\"\n", p->program);
      if (version != 0) 
        printfone("    version:      \"%s\"\n", p->version);
  }*/
}
#endif
