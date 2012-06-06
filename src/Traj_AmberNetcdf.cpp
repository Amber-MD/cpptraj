#ifdef BINTRAJ
// This file contains a collection of routines designed for reading
// netcdf trajectory files used with amber.
// Dan Roe 10-2008
// Original implementation of netcdf in Amber by Jon Mongan.
#include "netcdf.h"
#include "Traj_AmberNetcdf.h"
#include "CpptrajStdio.h"

// DEFINES 
#define NCREMD_DIMENSION "remd_dimension"
#define NCREMD_GROUPNUM "remd_groupnum"
#define NCREMD_DIMTYPE "remd_dimtype"
#define NCREMD_INDICES "remd_indices"

// CONSTRUCTOR
AmberNetcdf::AmberNetcdf() :
  Coord_(0),
  remd_dimension_(0),
  dimensionDID_(-1),
  groupnumVID_(-1),
  dimtypeVID_(-1),
  indicesVID_(-1),
  remd_groupnum_(0),
  remd_dimtype_(0),
  remd_indices_(0)
{
  //fprintf(stderr,"Amber Netcdf Constructor\n");
  // Netcdf files are always seekable
  seekable_=true;
}

// DESTRUCTOR
AmberNetcdf::~AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  this->closeTraj();
  if (Coord_!=0) delete[] Coord_;
  if (remd_groupnum_!=0) delete[] remd_groupnum_;
  if (remd_dimtype_!=0) delete[] remd_dimtype_;
  if (remd_indices_!=0) delete[] remd_indices_;
  // NOTE: Need to close file?
}

bool AmberNetcdf::ID_TrajFormat() {
  if ( GetNetcdfConventions( Name() ) == NC_AMBERTRAJ ) return true;
  return false;
} 

// AmberNetcdf::close()
/** Close netcdf file. Set ncid to -1 since it can change between open
  * and close calls.
  */
void AmberNetcdf::closeTraj() {
  NC_close();
}

// AmberNetcdf::openTraj()
/** Open up Netcdf file and set ncid. Variable and Dimension IDs are set up
  * by SetupRead / SetupWrite and will not change for a given file between
  * open and close calls.
  */
int AmberNetcdf::openTraj() {
  //fprintf(stdout,"DEBUG: AmberNetcdf::openTraj() called for %s, ncid=%i\n",BaseName(),ncid);
  // If already open, return
  if (Ncid()!=-1) return 0;

  switch (access_) {
    case READ :
      if ( NC_openRead( Name() ) != 0 ) {
        mprinterr("Error: Opening Netcdf file %s for reading.\n",BaseName()); 
        return 1;
      }
      break;
    
    case APPEND: 
    case WRITE:
      if ( NC_openWrite( Name() ) != 0 ) {
        mprinterr("Error: Opening Netcdf file %s for Write.\n",BaseName());
        return 1;
      }
  }

  if (debug_>0) rprintf("Successfully opened %s, ncid=%i\n",BaseName(),Ncid());
  if (debug_>1) NetcdfDebug();

  return 0;
}

// AmberNetcdf::setupTrajin()
/** Open the netcdf file, read all dimension and variable IDs, close.
  * Return the number of frames in the file. 
  */
int AmberNetcdf::setupTrajin(Topology* trajParm) {
  if (openTraj()) return -1;

  // Sanity check - Make sure this is a Netcdf trajectory
  if ( GetNetcdfConventions() != NC_AMBERTRAJ ) {
    mprinterr("Error: Netcdf file %s conventions do not include \"AMBER\"\n",BaseName());
    return -1;
  }
  // Get global attributes
  std::string attrText = GetAttrText("ConventionVersion");
  if ( attrText != "1.0") 
    mprintf("Warning: Netcdf file %s has ConventionVersion that is not 1.0 (%s)\n",
            BaseName(), attrText.c_str());
  // Get title
  if (title_.empty()) 
    title_ = GetAttrText("title");

  // Get Frame info
  if ( SetupFrame()!=0 ) return -1;

  // Setup Coordinates
  if ( SetupCoordinates()!=0 ) return -1;
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in NetCDF file %s (%i) does not\n",
              BaseName(),Ncatom());
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->Natom());
    return -1;
  }

  // Setup Time
  if ( SetupTime()!=0 ) return -1;

  // Box info
  // NOTE: If no box info found in parm should really try to determine correct
  //       box type from angles.
  int boxerr = SetupBox(boxAngle_, boxLength_);
  if (boxerr == 1)
    return 1;
  else if (boxerr == 0)
    hasBox_ = true;
  else
    hasBox_ = false;
 
  // Replica Temperatures - Allowed to fail gracefully
  if (SetupTemperature() == 0)
    hasTemperature_ = true;

/*  // Multi-d REMD info
  if ( nc_inq_dimid(ncid, NCREMD_DIMENSION, &dimensionDID) == NC_NOERR) {
    // Although this is a second call to dimid, makes for easier code
    if ( (dimensionDID = GetDimInfo(ncid, NCREMD_DIMENSION, &remd_dimension))==-1 )
      return -1;
    mprintf("    Netcdf file has multi-D REMD info, %i dimensions.\n",remd_dimension);
    // Ensure valid # dimensions
    if (remd_dimension < 1) {
      mprinterr("Error: Number of REMD dimensions is less than 1!\n");
      return -1;
    }
    // Start and count for groupnum and dimtype, allocate mem
    start[0]=0; start[1]=0; start[2]=0;
    count[0]=remd_dimension; count[1]=0; count[2]=0;
    remd_groupnum = new int[ remd_dimension ];
    remd_dimtype = new int[ remd_dimension ];
    remd_indices = new int[ remd_dimension ];
    // Get number of groups in each dimension
    if ( checkNCerr(nc_inq_varid(ncid, NCREMD_GROUPNUM, &groupnumVID),
                    "Getting group variable ID for each dimension.")!=0 ) return -1; 
    if ( checkNCerr(nc_get_vara_int(ncid, groupnumVID, start, count, remd_groupnum),
                    "Getting group numbers in each dimension.")!=0 ) return -1;
    // Get dimension types
    if ( checkNCerr(nc_inq_varid(ncid, NCREMD_DIMTYPE, &dimtypeVID),
                    "Getting dimension type variable ID for each dimension.")!=0 ) return -1;
    if ( checkNCerr(nc_get_vara_int(ncid, dimtypeVID, start, count, remd_dimtype),
                    "Getting dimension type in each dimension.")!=0 ) return -1;
    // Get VID for replica indices
    if ( checkNCerr(nc_inq_varid(ncid, NCREMD_INDICES, &indicesVID),
                    "Getting replica indices variable ID.")!=0 ) return -1;
    // Print info for each dimension
    for (int dim = 0; dim < remd_dimension; dim++)
      mprintf("\tDim %i: type %i (%i)\n",dim+1, remd_dimtype[dim], remd_groupnum[dim]);
  }*/

  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;


  // Amber Netcdf coords are float. Allocate a float array for converting
  // float to/from double.
  if (Coord_ != 0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  closeTraj();

  seekable_ = true;

  return Ncframe();
}

// AmberNetcdf::processWriteArgs()
int AmberNetcdf::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("remdtraj")) this->SetTemperature();
  return 0;
}

// AmberNetcdf::setupTrajout()
/** Create Netcdf file specified by filename and set up dimension and
  * variable IDs. 
  */
int AmberNetcdf::setupTrajout(Topology *trajParm, int NframesToWrite) {
  // Set up title
  if (title_.empty())
    title_.assign("Cpptraj Generated trajectory");

  if ( NC_create( Name(), NC_AMBERTRAJ, trajParm->Natom(), hasVelocity_,
                  hasBox_, hasTemperature_, true, title_ ) )
    return 1; 

  // Allocate memory and close the file. It will be reopened WRITE
  if (Coord_!=0) delete[] Coord_;
  Coord_ = new float[ Ncatom3() ];
  closeTraj();
  
  return 0;
}

// AmberNetcdf::readFrame()
/** Get the specified frame from amber netcdf file
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  */
int AmberNetcdf::readFrame(int set,double *X, double *V,double *box, double *T) {
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

  // Get replica indices
  if (indicesVID_!=-1) {
    start_[0] = set;
    start_[1] = 0;
    count_[0] = 1;
    count_[1] = remd_dimension_;
    if ( checkNCerr(nc_get_vara_int(ncid_, indicesVID_, start_, count_, remd_indices_)) ) {
      mprinterr("Error: Getting replica indices.\n");
      return 1;
    }
    mprintf("DEBUG:\tReplica indices:");
    for (int dim=0; dim < remd_dimension_; dim++) mprintf(" %i",remd_indices_[dim]);
    mprintf("\n");
  }

  // Read Coords 
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  if ( checkNCerr(nc_get_vara_float(ncid_, coordVID_, start_, count_, Coord_)) ) {
    mprinterr("Error: Getting frame %i\n",set);
    return 1;
  }
  // Read box info 
  if (hasBox_) {
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

  FloatToDouble(X, Coord_);

  return 0;
}

// AmberNetcdf::writeFrame() 
int AmberNetcdf::writeFrame(int set, double *X, double *V, double *box, double T) {

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
  if (hasBox_) {
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

// AmberNetcdf::info()
void AmberNetcdf::info() {
  mprintf("is a NetCDF AMBER trajectory"
            //(p->isVelocity ? " and velocities" : "")
         );

  if (hasTemperature_) mprintf(" with replica temperatures");

  /*if (debug_ > 2) {
      if (title != NULL)
        printfone("    title:        \"%s\"\n", p->title);
      if (application != NULL)  
        printfone("    application:  \"%s\"\n", p->application);
      if (program != NULL) 
        printfone("    program:      \"%s\"\n", p->program);
      if (version != NULL) 
        printfone("    version:      \"%s\"\n", p->version);
  }*/
}
#endif
