#ifdef BINTRAJ
#  include <netcdf.h>
#endif
#include "NetcdfFile.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "Version.h"

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(const char* fname) {
  NCTYPE nctype = NC_UNKNOWN;
#ifdef BINTRAJ
  // NOTE: Do not use checkNCerr so this fails silently. Allows routine to
  //       be used in file autodetection.
  if ( nc_open( fname, NC_NOWRITE, &ncid_ ) != NC_NOERR )
    return NC_UNKNOWN;
  nctype = GetNetcdfConventions();
  NC_close();
#else
  mprintf("Error: Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#endif
  return nctype;
}

#ifdef BINTRAJ
// DEFINES
#define NCENSEMBLE "ensemble"
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_LENGTHS "cell_lengths"
#define NCCELL_ANGULAR "cell_angular"
#define NCCELL_ANGLES "cell_angles"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCFRC "forces"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5
#define NCREMD_DIMENSION "remd_dimension"
#define NCREMD_DIMTYPE "remd_dimtype"
#define NCREMD_INDICES "remd_indices"
#define NCREMD_REPIDX "remd_repidx"
#define NCREMD_CRDIDX "remd_crdidx"
#define NCEPTOT "eptot"
#define NCBINS "bins"

// CONSTRUCTOR
NetcdfFile::NetcdfFile() :
  ncid_(-1),
  ncframe_(-1),
  TempVID_(-1),
  coordVID_(-1),
  velocityVID_(-1),
  frcVID_(-1),
  cellAngleVID_(-1),
  cellLengthVID_(-1),
  timeVID_(-1),
  remd_dimension_(0),
  indicesVID_(-1),
  repidxVID_(-1),
  crdidxVID_(-1),
  ncdebug_(0),
  frameDID_(-1),
  atomDID_(-1),
  ncatom_(-1),
  ncatom3_(-1),
  spatialDID_(-1),
  labelDID_(-1),
  cell_spatialDID_(-1),
  cell_angularDID_(-1),
  spatialVID_(-1),
  cell_spatialVID_(-1),
  cell_angularVID_(-1)
{
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 0;
  count_[1] = 0;
  count_[2] = 0;
}

// NetcdfFile::GetAttrText()
/** Get the information about a netcdf attribute with given vid and 
  * attribute text.
  * Since there is no guarantee that null char at the end of retrieved string
  * append one.
  */
std::string NetcdfFile::GetAttrText(int vid, const char *attribute) {
  size_t attlen;
  std::string attrOut;
  // Get attr length
  if ( checkNCerr(nc_inq_attlen(ncid_, vid, attribute, &attlen)) ) {
    mprintf("Warning: Getting length for attribute '%s'\n",attribute);
    return attrOut;
  }
  // Allocate space for attr text, plus one for null char
  char *attrText = new char[ (attlen + 1) ];
  // Get attr text
  if ( checkNCerr(nc_get_att_text(ncid_, vid, attribute, attrText)) ) {
    mprintf("Warning: Getting attribute text for '%s'\n",attribute);
    delete[] attrText;
    return attrOut;
  }
  // Append null char - NECESSARY?
  attrText[attlen]='\0';
  attrOut.assign( attrText );
  delete[] attrText;

  return attrOut;
}

// NetcdfFile::GetAttrText()
/** Get information about a netcdf global attribute. */
std::string NetcdfFile::GetAttrText(const char *attribute) {
  return GetAttrText(NC_GLOBAL, attribute);
}

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions() {
  NCTYPE nctype = NC_UNKNOWN;
  std::string attrText = GetAttrText(NC_GLOBAL, "Conventions");
  if (attrText == "AMBERENSEMBLE")
    nctype = NC_AMBERENSEMBLE;
  else if (attrText == "AMBER")
    nctype = NC_AMBERTRAJ;
  else if (attrText == "AMBERRESTART")
    nctype = NC_AMBERRESTART;
  else if (attrText.empty()) 
    mprinterr("Error: Could not get conventions from Netcdf file.\n");
  else {
    mprinterr("Error: Netcdf file: Unrecognized conventions \"%s\".\n",
              attrText.c_str());
    mprinterr("Error:   Expected \"AMBER\", \"AMBERRESTART\", or \"AMBERENSEMBLE\".\n");
  }
  return nctype;
}

// NetcdfFile::GetDimInfo()
/** Return the dimension ID of a given attribute in netcdf file ncid.
  * Also set dimension length.
  */
int NetcdfFile::GetDimInfo(const char *attribute, int *length) {
  int dimID;
  size_t slength = 0;

  *length = 0;
  // Get dimid 
  if ( checkNCerr(nc_inq_dimid(ncid_, attribute, &dimID)) ) {
    mprinterr("Error: Getting dimID for attribute %s\n",attribute);
    return -1;
  }
  // get Dim length 
  if ( checkNCerr(nc_inq_dimlen(ncid_, dimID, &slength)) ) {
    mprinterr("Error: Getting length for attribute %s\n",attribute);
    return -1;
  }
  *length = (int) slength;
  return dimID;
}

// NetcdfFile::SetupFrameDim()
/** Get the frame dimension ID and # of frames (ncframe). */
int NetcdfFile::SetupFrameDim() {
  frameDID_ = GetDimInfo( NCFRAME, &ncframe_ );
  if (frameDID_==-1) return 1;
  return 0;
}

/** Get the ensemble dimension ID and size. */
int NetcdfFile::SetupEnsembleDim() {
  int ensembleSize = 0;
  ensembleDID_ = GetDimInfo( NCENSEMBLE, &ensembleSize );
  if (ensembleDID_ == -1) return 0;
  return ensembleSize;
}

// NetcdfFile::SetupCoordsVelo()
/** Setup ncatom, ncatom3, atomDID, coordVID, spatialDID, spatialVID,
  * velocityVID, frcVID. Check units and spatial dimensions.
  */
int NetcdfFile::SetupCoordsVelo(bool useVelAsCoords) {
  int spatial;
  atomDID_ = GetDimInfo(NCATOM, &ncatom_);
  if (atomDID_==-1) return 1;
  ncatom3_ = ncatom_ * 3;
  // Get coord info
  coordVID_ = -1;
  if ( nc_inq_varid(ncid_, NCCOORDS, &coordVID_) == NC_NOERR ) {
    if (ncdebug_ > 0) mprintf("\tNetcdf file has coordinates.\n");
    std::string attrText = GetAttrText(coordVID_, "units");
    if (attrText!="angstrom")
      mprintf("Warning: Netcdf file has length units of %s - expected angstrom.\n",
              attrText.c_str());
  }
  // Get spatial info
  spatialDID_ = GetDimInfo(NCSPATIAL, &spatial);
  if (spatialDID_==-1) return 1;
  if (spatial!=3) {
    mprinterr("Error: Expected 3 spatial dimensions, got %i\n",spatial);
    return 1;
  }
  if ( checkNCerr(nc_inq_varid(ncid_, NCSPATIAL, &spatialVID_)) ) {
    mprintf("Warning: Could not get spatial VID. File may not be Amber NetCDF compliant.\n");
    mprintf("Warning: Assuming spatial variables are 'x', 'y', 'z'\n");
  } else {
    start_[0] = 0;
    count_[0] = 3;
    char xyz[3];
    if (checkNCerr(nc_get_vara_text(ncid_, spatialVID_, start_, count_, xyz))) {
      mprinterr("Error: Getting spatial variables.\n");
      return 1;
    }
    if (xyz[0] != 'x' || xyz[1] != 'y' || xyz[2] != 'z') {
      mprinterr("Error: NetCDF spatial variables are '%c', '%c', '%c', not 'x', 'y', 'z'\n",
                xyz[0], xyz[1], xyz[2]);
      return 1;
    }
  }
  // Get velocity info
  velocityVID_ = -1;
  if ( nc_inq_varid(ncid_, NCVELO, &velocityVID_) == NC_NOERR ) {
    if (ncdebug_>0) mprintf("\tNetcdf file has velocities.\n");
  }
  // Return a error if no coords and no velocity
  if ( coordVID_ == -1 && velocityVID_ == -1 ) {
    mprinterr("Error: NetCDF file has no coords and no velocities.\n");
    return 1;
  }
  // If using velocities as coordinates, swap them now.
  if (useVelAsCoords) {
    if (velocityVID_ == -1) {
      mprinterr("Error: Cannot use velocities as coordinates; no velocities present.\n");
      return 1;
    }
    mprintf("\tUsing velocities as coordinates.\n");
    coordVID_ = velocityVID_;
    velocityVID_ = -1;
  }
  // Get force info
  frcVID_ = -1;
  if ( nc_inq_varid(ncid_, NCFRC, &frcVID_) == NC_NOERR ) {
    if (ncdebug_>0) mprintf("\tNetcdf file has forces.\n");
  }
  // Get overall replica and coordinate indices
  crdidxVID_ = -1;
  if ( nc_inq_varid(ncid_, NCREMD_REPIDX, &repidxVID_) == NC_NOERR ) {
      if (ncdebug_>0) mprintf("\tNetcdf file has overall replica indices.\n");
    if ( checkNCerr(nc_inq_varid(ncid_, NCREMD_CRDIDX, &crdidxVID_)) ) {
      mprinterr("Error: Getting overall coordinate index variable ID.\n");
      return 1;
    }
  } else
    repidxVID_ = -1;
  return 0;
}

// NetcdfFile::SetupTime()
/** Determine if Netcdf file contains time; set up timeVID and check units. */
int NetcdfFile::SetupTime() {
  if ( nc_inq_varid(ncid_, NCTIME, &timeVID_) == NC_NOERR ) {
    std::string attrText = GetAttrText(timeVID_, "units");
    if (attrText!="picosecond")
      mprintf("Warning: NetCDF file has time units of %s - expected picosecond.\n",
              attrText.c_str());
    // Check for time values which have NOT been filled, which was possible
    // with netcdf trajectories created by older versions of ptraj/cpptraj.
    if (ncframe_ > 0 && GetNetcdfConventions() == NC_AMBERTRAJ) {
      float time;
      start_[0] = 0; count_[0] = 1;
      if (checkNCerr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time))) {
        mprinterr("Error: Getting time value for NetCDF file.\n");
        return -1;
      }
      if (time == NC_FILL_FLOAT) {
        mprintf("Warning: NetCDF file time variable defined but empty. Disabling.\n");
        timeVID_ = -1;
      }
    }
    return 0;
  }
  timeVID_=-1;
  return 1;
}

// NetcdfFile::SetupTemperature()
/** Determine if Netcdf file contains temperature; set up temperature VID. */
int NetcdfFile::SetupTemperature() {
  if ( nc_inq_varid(ncid_,NCTEMPERATURE,&TempVID_) == NC_NOERR ) {
    if (ncdebug_>0) mprintf("    Netcdf file has replica temperatures.\n");
    return 0;
  } 
  TempVID_=-1;
  return 1;
}

// NetcdfFile::SetupMultiD()
/** Determine if Netcdf file contains multi-D REMD info. If so set the
  * number of replica dimensions (remd_dimension_) and figure out
  * the dimension types (remdDim)
  */
int NetcdfFile::SetupMultiD(ReplicaDimArray& remdDim) {
  int dimensionDID;

  if ( nc_inq_dimid(ncid_, NCREMD_DIMENSION, &dimensionDID) != NC_NOERR)
    return 1;
 
  // Although this is a second call to dimid, makes for easier code
  if ( (dimensionDID = GetDimInfo(NCREMD_DIMENSION, &remd_dimension_))==-1 )
    return -1;
  if (ncdebug_ > 0)
    mprintf("\tNetcdf file has multi-D REMD info, %i dimensions.\n",remd_dimension_);
  // Ensure valid # dimensions
  if (remd_dimension_ < 1) {
    mprinterr("Error: Number of REMD dimensions is less than 1!\n");
    return -1;
  }
  // Start and count for groupnum and dimtype, allocate mem
  start_[0]=0; 
  start_[1]=0; 
  start_[2]=0;
  count_[0]=remd_dimension_; 
  count_[1]=0; 
  count_[2]=0;
  int* remd_dimtype = new int[ remd_dimension_ ];
  // Get dimension types
  int dimtypeVID;
  if ( checkNCerr(nc_inq_varid(ncid_, NCREMD_DIMTYPE, &dimtypeVID)) ) {
    mprinterr("Error: Getting dimension type variable ID for each dimension.\n");
    return -1;
  }
  if ( checkNCerr(nc_get_vara_int(ncid_, dimtypeVID, start_, count_, remd_dimtype)) ) {
    mprinterr("Error: Getting dimension type in each dimension.\n");
    return -1;
  }
  // Get VID for replica indices
  if ( checkNCerr(nc_inq_varid(ncid_, NCREMD_INDICES, &indicesVID_)) ) {
    mprinterr("Error: Getting replica indices variable ID.\n");
    return -1;
  }
  // Print info for each dimension
  for (int dim = 0; dim < remd_dimension_; ++dim)
    remdDim.AddRemdDimension( remd_dimtype[dim] );
  delete[] remd_dimtype;
  return 0; 
}

// NetcdfFile::SetupBox()
/** \return 0 on success, 1 on error, -1 for no box coords. */
// TODO: Use Box class
int NetcdfFile::SetupBox(double* boxIn, NCTYPE typeIn) {
  boxIn[0] = 0.0;
  boxIn[1] = 0.0;
  boxIn[2] = 0.0;
  boxIn[3] = 0.0;
  boxIn[4] = 0.0;
  boxIn[5] = 0.0;
  if ( nc_inq_varid(ncid_,NCCELL_LENGTHS,&cellLengthVID_)==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(ncid_,NCCELL_ANGLES,&cellAngleVID_)) ) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (ncdebug_>0) mprintf("  Netcdf Box information found.\n");
    // If present, get box angles. These will be used to determine the box 
    // type in TrajectoryFile.
    start_[0]=0; 
    start_[1]=0; 
    start_[2]=0;
    start_[3]=0;
    switch (typeIn) {
      case NC_AMBERRESTART:
        count_[0]=3;
        count_[1]=0;
        count_[2]=0;
        break;
      case NC_AMBERTRAJ:
        count_[0]=1; 
        count_[1]=3;
        count_[2]=0;
        break;
      case NC_AMBERENSEMBLE:
        count_[0]=1; // NOTE: All ensemble members must have same box type
        count_[1]=1; // TODO: Check all members?
        count_[2]=3;
        break;
      case NC_UNKNOWN: return 1; // Sanity check
    }
    count_[3]=0;
    if ( checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, boxIn )) )
    {
      mprinterr("Error: Getting cell lengths.\n");
      return 1;
    }
    if ( checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, boxIn+3)) )
    {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (ncdebug_ > 0) mprintf("\tNetcdf Box: XYZ={%f %f %f} ABG={%f %f %f}\n",
                              boxIn[0], boxIn[1], boxIn[2], boxIn[3], boxIn[4], boxIn[5]);
    return 0;
  }
  // No box information
  return -1;
}

// NetcdfFile::checkNCerr()
bool NetcdfFile::checkNCerr(int ncerr) {
  if ( ncerr != NC_NOERR ) {
    mprintf("%s\n", nc_strerror(ncerr));
    return true;
  }
  return false;
}

// NetcdfFile::NetcdfDebug()
void NetcdfFile::NetcdfDebug() {
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  char varname[128];
  // ncid:    NetCDF ID, from a previous call to nc open or nc create.
  // ndimsp:  Pointer to location for returned number of dimensions defined for 
  //         this netCDF dataset.
  // nvarsp:  Pointer to location for returned number of variables defined for 
  //         this netCDF dataset.
  // ngattsp: Pointer to location for returned number of global attributes 
  //         defined for this netCDF dataset.
  // unlimdimidp: 
  //  Pointer to location for returned ID of the unlimited dimension, if 
  //  there is one for this netCDF dataset. If no unlimited length 
  //  dimension has been defined, -1 is returned.
  mprintf("========== BEG. NETCDF DEBUG ==========\n");
  int err = nc_inq(ncid_, &ndimsp, &nvarsp, &ngattsp, &unlimdimidp);
  mprintf("nc_inq returned %i\n",err);
  if (err==NC_NOERR)
    mprintf("ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            ndimsp,nvarsp,ngattsp,unlimdimidp);
  else
    mprintf("NETCDF Error occurred.\n");
  // Print name of each variable defined in netcdf file
  mprintf("NC VARIABLES:\n");
  for (int i=0; i<nvarsp; i++) {
    err=nc_inq_varname(ncid_,i,varname);
    mprintf("  Var %i - ",i);
    if (err==NC_NOERR)
      mprintf("%s\n",varname);
    else
      mprintf("NETCDF Error occured.\n");
  }
  mprintf("==========  END NETCDF DEBUG ==========\n");
}

// NetcdfFile::NC_openRead()
int NetcdfFile::NC_openRead(std::string const& Name) {
  if (Name.empty()) return 1;
  if ( checkNCerr( nc_open( Name.c_str(), NC_NOWRITE, &ncid_ ) ) )
    return 1;
  return 0;
}

// NetcdfFile::NC_close()
void NetcdfFile::NC_close() {
  if (ncid_ == -1) return;
  bool err = checkNCerr( nc_close(ncid_) );
  if (ncdebug_ > 0 && !err)
    mprintf("Successfully closed ncid %i\n",ncid_);
  ncid_ = -1;
}

// NetcdfFile::NC_openWrite()
int NetcdfFile::NC_openWrite(std::string const& Name) {
  if (Name.empty()) return 1;
  if ( checkNCerr( nc_open( Name.c_str(), NC_WRITE, &ncid_ ) ) )
    return 1;
  return 0;
}

// NetcdfFile::NC_defineTemperature()
int NetcdfFile::NC_defineTemperature(int* dimensionID, int NDIM) {
  if (checkNCerr(nc_def_var(ncid_,NCTEMPERATURE,NC_DOUBLE,NDIM,dimensionID,&TempVID_))) {
    mprinterr("NetCDF error on defining temperature.\n");
    return 1;
  }
  if (checkNCerr(nc_put_att_text(ncid_,TempVID_,"units",6,"kelvin"))) {
    mprinterr("NetCDF error on defining temperature units.\n");
    return 1;
  }
  return 0;
}

// NetcdfFile::NC_createReservoir()
int NetcdfFile::NC_createReservoir(bool hasBins, double reservoirT, int iseed,
                                   int& eptotVID, int& binsVID) 
{
  int dimensionID[1];
  dimensionID[0] = frameDID_;
  if (ncid_ == -1 || dimensionID[0] == -1) return 1;
  // Place file back in define mode
  if ( checkNCerr( nc_redef( ncid_ ) ) ) return 1;
  // Define eptot, bins, temp0
  if ( checkNCerr( nc_def_var(ncid_, NCEPTOT, NC_DOUBLE, 1, dimensionID, &eptotVID)) ) {
    mprinterr("Error: defining eptot variable ID.\n");
    return 1;
  }
  if (hasBins) {
    if ( checkNCerr( nc_def_var(ncid_, NCBINS, NC_INT, 1, dimensionID, &binsVID)) ) {
      mprinterr("Error: defining bins variable ID.\n");
      return 1;
    }
  } else
    binsVID = -1;
  if (NC_defineTemperature(dimensionID, 0)) return 1;
  // Random seed, make global
  if ( checkNCerr( nc_put_att_int(ncid_, NC_GLOBAL, "iseed", NC_INT, 1, &iseed) ) ) {
    mprinterr("Error: setting random seed attribute.\n");
    return 1;
  }
  // End definitions
  if (checkNCerr(nc_enddef(ncid_))) {
    mprinterr("NetCDF error on ending definitions.");
    return 1;
  }
  // Write temperature
  if (checkNCerr(nc_put_var_double(ncid_,TempVID_,&reservoirT)) ) {
    mprinterr("Error: Writing reservoir temperature.\n");
    return 1;
  }
  return 0;
}

// NetcdfFile::NC_create()
int NetcdfFile::NC_create(std::string const& Name, NCTYPE type, int natomIn,
                          CoordinateInfo const& coordInfo, std::string const& title) 
{
  if (Name.empty()) return 1;
  int dimensionID[NC_MAX_VAR_DIMS];
  int NDIM;
  nc_type dataType;

  if (ncdebug_>1)
    mprintf("DEBUG: NC_create: %s  natom=%i V=%i F=%i box=%i  temp=%i  time=%i\n",
            Name.c_str(),natomIn,(int)coordInfo.HasVel(),
            (int)coordInfo.HasForce(),(int)coordInfo.HasBox(),
            (int)coordInfo.HasTemp(),(int)coordInfo.HasTime());

  if ( checkNCerr( nc_create( Name.c_str(), NC_64BIT_OFFSET, &ncid_) ) )
    return 1;

  ncatom_ = natomIn;
  ncatom3_ = ncatom_ * 3;
  
  // Set number of dimensions based on file type
  switch (type) {
    case NC_AMBERENSEMBLE:
      NDIM = 4;
      dataType = NC_FLOAT;
      break;
    case NC_AMBERTRAJ: 
      NDIM = 3;
      dataType = NC_FLOAT;
      break;
    case NC_AMBERRESTART: 
      NDIM = 2; 
      dataType = NC_DOUBLE;
      break;
    default:
      mprinterr("Error: NC_create (%s): Unrecognized type (%i)\n",Name.c_str(),(int)type);
      return 1;
  }

  if (type == NC_AMBERENSEMBLE) {
    // Ensemble dimension for ensemble
    if (coordInfo.EnsembleSize() < 1) {
      mprinterr("Internal Error: NetcdfFile: ensembleSize < 1\n");
      return 1;
    }
    if ( checkNCerr(nc_def_dim(ncid_, NCENSEMBLE, coordInfo.EnsembleSize(), &ensembleDID_)) ) {
      mprinterr("Error: Defining ensemble dimension.\n");
      return 1;
    }
    dimensionID[1] = ensembleDID_;
  }
  ncframe_ = 0;
  if (type == NC_AMBERTRAJ || type == NC_AMBERENSEMBLE) {
    // Frame dimension for traj
    if ( checkNCerr( nc_def_dim( ncid_, NCFRAME, NC_UNLIMITED, &frameDID_)) ) {
      mprinterr("Error: Defining frame dimension.\n");
      return 1;
    }
    // Since frame is UNLIMITED, it must be lowest dim.
    dimensionID[0] = frameDID_;
  }
  // Time variable and units
  if (coordInfo.HasTime()) {
    if ( checkNCerr( nc_def_var(ncid_, NCTIME, dataType, NDIM-2, dimensionID, &timeVID_)) ) {
      mprinterr("Error: Defining time variable.\n");
      return 1;
    }
    if ( checkNCerr( nc_put_att_text(ncid_, timeVID_, "units", 10, "picosecond")) ) {
      mprinterr("Error: Writing time VID units.\n");
      return 1;
    }
  }
  // Spatial dimension and variable
  if ( checkNCerr( nc_def_dim( ncid_, NCSPATIAL, 3, &spatialDID_)) ) {
    mprinterr("Error: Defining spatial dimension.\n");
    return 1;
  }
  dimensionID[0] = spatialDID_;
  if ( checkNCerr( nc_def_var( ncid_, NCSPATIAL, NC_CHAR, 1, dimensionID, &spatialVID_)) ) {
    mprinterr("Error: Defining spatial variable.\n"); 
    return 1;
  }
  // Atom dimension
  if ( checkNCerr( nc_def_dim( ncid_, NCATOM, ncatom_, &atomDID_)) ) {
    mprinterr("Error: Defining atom dimension.\n");
    return 1;
  }
  // Setup dimensions for Coords/Velocity
  // NOTE: THIS MUST BE MODIFIED IF NEW TYPES ADDED
  if (type == NC_AMBERENSEMBLE) {
    dimensionID[0] = frameDID_;
    dimensionID[1] = ensembleDID_;
    dimensionID[2] = atomDID_;
    dimensionID[3] = spatialDID_;
  } else if (type == NC_AMBERTRAJ) {
    dimensionID[0] = frameDID_;
    dimensionID[1] = atomDID_;
    dimensionID[2] = spatialDID_;
  } else {
    dimensionID[0] = atomDID_;
    dimensionID[1] = spatialDID_;
  }
  // Coord variable
  if ( checkNCerr( nc_def_var( ncid_, NCCOORDS, dataType, NDIM, dimensionID, &coordVID_)) ) {
    mprinterr("Error: Defining coordinates variable.\n");
    return 1;
  }
  if ( checkNCerr( nc_put_att_text( ncid_, coordVID_, "units", 8, "angstrom")) ) {
    mprinterr("Error: Writing coordinates variable units.\n");
    return 1;
  }
  // Velocity variable
  if (coordInfo.HasVel()) {
    if ( checkNCerr( nc_def_var( ncid_, NCVELO, dataType, NDIM, dimensionID, &velocityVID_)) ) {
      mprinterr("Error: Defining velocities variable.\n");
      return 1;
    }
    if ( checkNCerr( nc_put_att_text( ncid_, velocityVID_, "units", 19, "angstrom/picosecond")) )
    {
      mprinterr("Error: Writing velocities variable units.\n");
      return 1;
    }
    if ( checkNCerr( nc_put_att_double( ncid_, velocityVID_, "scale_factor", NC_DOUBLE, 1, 
                                        &Constants::AMBERTIME_TO_PS)) )
    {
      mprinterr("Error: Writing velocities scale factor.\n");
      return 1;
    }
  }
  // Force variable
  if (coordInfo.HasForce()) {
    if ( checkNCerr( nc_def_var( ncid_, NCFRC, dataType, NDIM, dimensionID, &frcVID_)) ) {
      mprinterr("Error: Defining forces variable\n");
      return 1;
    }
    if ( checkNCerr( nc_put_att_text( ncid_, frcVID_, "units", 25, "kilocalorie/mole/angstrom")) )
    {
      mprinterr("Error: Writing forces variable units.\n");
      return 1;
    }
  }
  // Replica Temperature
  if (coordInfo.HasTemp()) {
    // NOTE: Setting dimensionID should be OK for Restart, will not be used.
    dimensionID[0] = frameDID_;
    if ( NC_defineTemperature( dimensionID, NDIM-2 ) ) return 1;
  }
  // Replica indices
  int remDimTypeVID = -1;
  if (coordInfo.HasReplicaDims()) {
    // Define number of replica dimensions
    remd_dimension_ = coordInfo.ReplicaDimensions().Ndims();
    int remDimDID = -1;
    if ( checkNCerr(nc_def_dim(ncid_, NCREMD_DIMENSION, remd_dimension_, &remDimDID)) ) {
      mprinterr("Error: Defining replica indices dimension.\n");
      return 1;
    }
    dimensionID[0] = remDimDID;
    // For each dimension, store the type
    if ( checkNCerr(nc_def_var(ncid_, NCREMD_DIMTYPE, NC_INT, 1, dimensionID, &remDimTypeVID)) ) 
    {
      mprinterr("Error: Defining replica dimension type variable.\n");
      return 1;
    }
    // Need to store the indices of replica in each dimension each frame
    // NOTE: THIS MUST BE MODIFIED IF NEW TYPES ADDED
    if (type == NC_AMBERENSEMBLE) {
      dimensionID[0] = frameDID_;
      dimensionID[1] = ensembleDID_;
      dimensionID[2] = remDimDID;
    } else if (type == NC_AMBERTRAJ) {
      dimensionID[0] = frameDID_;
      dimensionID[1] = remDimDID;
    } else {
      dimensionID[0] = remDimDID;
    }
    if (checkNCerr(nc_def_var(ncid_, NCREMD_INDICES, NC_INT, NDIM-1, dimensionID, &indicesVID_)))
    {
      mprinterr("Error: Defining replica indices variable ID.\n");
      return 1;
    }
    // TODO: Determine if groups are really necessary for restarts. If not, 
    // remove from AmberNetcdf.F90.
  }
  // Box Info
  if (coordInfo.HasBox()) {
    // Cell Spatial
    if ( checkNCerr( nc_def_dim( ncid_, NCCELL_SPATIAL, 3, &cell_spatialDID_)) ) {
      mprinterr("Error: Defining cell spatial dimension.\n");
      return 1;
    }
    dimensionID[0] = cell_spatialDID_;
    if ( checkNCerr( nc_def_var(ncid_, NCCELL_SPATIAL, NC_CHAR, 1, dimensionID, &cell_spatialVID_)))
    {
      mprinterr("Error: Defining cell spatial variable.\n");
      return 1;
    }
    // Cell angular
    if ( checkNCerr( nc_def_dim( ncid_, NCLABEL, NCLABELLEN, &labelDID_)) ) {
      mprinterr("Error: Defining label dimension.\n");
      return 1;
    }
    if ( checkNCerr( nc_def_dim( ncid_, NCCELL_ANGULAR, 3, &cell_angularDID_)) ) {
      mprinterr("Error: Defining cell angular dimension.\n"); 
      return 1;
    }
    dimensionID[0] = cell_angularDID_;
    dimensionID[1] = labelDID_;
    if ( checkNCerr( nc_def_var( ncid_, NCCELL_ANGULAR, NC_CHAR, 2, dimensionID, 
                                 &cell_angularVID_)) )
    {
      mprinterr("Error: Defining cell angular variable.\n");
      return 1;
    }
    // Setup dimensions for Box
    // NOTE: This must be modified if more types added
    int boxdim;
    if (type == NC_AMBERENSEMBLE) {
      dimensionID[0] = frameDID_;
      dimensionID[1] = ensembleDID_;
      boxdim = 2;
    } else if (type == NC_AMBERTRAJ) {
      dimensionID[0] = frameDID_;
      boxdim = 1;
    } else {
      boxdim = 0;
    }
    dimensionID[boxdim] = cell_spatialDID_;
    if ( checkNCerr( nc_def_var( ncid_, NCCELL_LENGTHS, NC_DOUBLE, NDIM-1, dimensionID,
                                 &cellLengthVID_)) )
    {
      mprinterr("Error: Defining cell length variable.\n"); 
      return 1;
    }
    if ( checkNCerr( nc_put_att_text( ncid_, cellLengthVID_, "units", 8, "angstrom")) ) {
      mprinterr("Error: Writing cell length variable units.\n");
      return 1;
    }
    dimensionID[boxdim] = cell_angularDID_;
    if ( checkNCerr( nc_def_var( ncid_, NCCELL_ANGLES, NC_DOUBLE, NDIM-1, dimensionID,
                                 &cellAngleVID_)) )
    {
      mprinterr("Error: Defining cell angle variable.\n");
      return 1;
    }
    if ( checkNCerr( nc_put_att_text( ncid_, cellAngleVID_, "units", 6, "degree")) ) {
      mprinterr("Error: Writing cell angle variable units.\n");
      return 1;
    }
  }

  // Attributes
  if (checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"title",title.size(),title.c_str())) ) {
    mprinterr("Error: Writing title.\n");
    return 1;
  }
  if (checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"application",5,"AMBER")) ) {
    mprinterr("Error: Writing application.\n");
    return 1;
  }
  if (checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"program",7,"cpptraj")) ) {
    mprinterr("Error: Writing program.\n");
    return 1;
  }
  if (checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"programVersion",
                                 NETCDF_VERSION_STRLEN, NETCDF_VERSION_STRING)) ) 
  {
    mprinterr("Error: Writing program version.\n");
    return 1;
  }
  // TODO: Make conventions a static string
  bool errOccurred = false;
  if ( type == NC_AMBERENSEMBLE )
    errOccurred = checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"Conventions",13,"AMBERENSEMBLE"));
  else if ( type == NC_AMBERTRAJ )
    errOccurred = checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"Conventions",5,"AMBER"));
  else
    errOccurred = checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"Conventions",12,"AMBERRESTART"));
  if (errOccurred) {
    mprinterr("Error: Writing conventions.\n");
    return 1;
  }
  if (checkNCerr(nc_put_att_text(ncid_,NC_GLOBAL,"ConventionVersion",3,"1.0")) ) {
    mprinterr("Error: Writing conventions version.\n");
    return 1;
  }
  
  // Set fill mode
  if (checkNCerr(nc_set_fill(ncid_, NC_NOFILL, dimensionID))) {
    mprinterr("Error: NetCDF setting fill value.\n");
    return 1;
  }

  // End netcdf definitions
  if (checkNCerr(nc_enddef(ncid_))) {
    mprinterr("NetCDF error on ending definitions.");
    return 1;
  }

  // Specify spatial dimension labels
  start_[0] = 0;
  count_[0] = 3;
  char xyz[3];
  xyz[0] = 'x'; 
  xyz[1] = 'y'; 
  xyz[2] = 'z';
  if (checkNCerr(nc_put_vara_text(ncid_, spatialVID_, start_, count_, xyz))) {
    mprinterr("Error on NetCDF output of spatial VID 'x', 'y' and 'z'");
    return 1;
  }
  if ( coordInfo.HasBox() ) {
    xyz[0] = 'a'; 
    xyz[1] = 'b'; 
    xyz[2] = 'c';
    if (checkNCerr(nc_put_vara_text(ncid_, cell_spatialVID_, start_, count_, xyz))) {
      mprinterr("Error on NetCDF output of cell spatial VID 'a', 'b' and 'c'");
      return 1;
    }
    char abc[15] = { 'a', 'l', 'p', 'h', 'a',
                     'b', 'e', 't', 'a', ' ',
                     'g', 'a', 'm', 'm', 'a' };
    start_[0] = 0; 
    start_[1] = 0;
    count_[0] = 3; 
    count_[1] = NCLABELLEN;
    if (checkNCerr(nc_put_vara_text(ncid_, cell_angularVID_, start_, count_, abc))) {
      mprinterr("Error on NetCDF output of cell angular VID 'alpha', 'beta ' and 'gamma'");
      return 1;
    }
  }

  // Store the type of each replica dimension.
  if (coordInfo.HasReplicaDims()) {
    ReplicaDimArray const& remdDim = coordInfo.ReplicaDimensions();
    start_[0] = 0;
    count_[0] = remd_dimension_;
    int* tempDims = new int[ remd_dimension_ ];
    for (int i = 0; i < remd_dimension_; ++i)
      tempDims[i] = remdDim[i];
    if (checkNCerr(nc_put_vara_int(ncid_, remDimTypeVID, start_, count_, tempDims))) {
      mprinterr("Error: writing replica dimension types.\n");
      delete[] tempDims;
      return 1;
    }
    delete[] tempDims;
  }

  return 0;
}

// FloatToDouble()
/** Convert float coords to double coords
  * NOTE: natom3 needs to match up with size of Coord!
  */
void NetcdfFile::FloatToDouble(double* X, const float* Coord) {
  for (int i=0; i<ncatom3_; ++i)
    X[i]=(double) Coord[i];
}

// DoubleToFloat()
/** Convert double coords to float coords
  */
void NetcdfFile::DoubleToFloat(float* Coord, const double* X) {
  for (int i=0; i<ncatom3_; ++i)
    Coord[i]=(float) X[i];
}

void NetcdfFile::WriteIndices() const {
  mprintf("DBG: Start={%zu, %zu, %zu, %zu} Count={%zu, %zu, %zu, %zu}\n",
         start_[0], start_[1], start_[2], start_[3],
         count_[0], count_[1], count_[2], count_[3]);
}

void NetcdfFile::WriteVIDs() const {
  rprintf("TempVID_=%i  coordVID_=%i  velocityVID_=%i frcVID_=%i  cellAngleVID_=%i"
          "  cellLengthVID_=%i  indicesVID_=%i\n",
          TempVID_, coordVID_, velocityVID_, frcVID_, cellAngleVID_, cellLengthVID_, indicesVID_);
}

#ifdef MPI
void NetcdfFile::Sync(Parallel::Comm const& commIn) {
  int nc_vars[23];
  if (commIn.Master()) {
    nc_vars[0] = ncframe_;
    nc_vars[1] = TempVID_;
    nc_vars[2] = coordVID_;
    nc_vars[3] = velocityVID_;
    nc_vars[4] = frcVID_;
    nc_vars[5] = cellAngleVID_;
    nc_vars[6] = cellLengthVID_;
    nc_vars[7] = timeVID_;
    nc_vars[8] = remd_dimension_;
    nc_vars[9] = indicesVID_;
    nc_vars[10] = ncdebug_;
    nc_vars[11] = ensembleDID_;
    nc_vars[12] = frameDID_;
    nc_vars[13] = atomDID_;
    nc_vars[14] = ncatom_;
    nc_vars[15] = ncatom3_;
    nc_vars[16] = spatialDID_;
    nc_vars[17] = labelDID_;
    nc_vars[18] = cell_spatialDID_;
    nc_vars[19] = cell_angularDID_;
    nc_vars[20] = spatialVID_;
    nc_vars[21] = cell_spatialVID_;
    nc_vars[22] = cell_angularVID_;
  }
  commIn.MasterBcast( nc_vars, 23, MPI_INT );
  if (!commIn.Master()) {
    ncframe_ = nc_vars[0];
    TempVID_ = nc_vars[1];
    coordVID_ = nc_vars[2];
    velocityVID_ = nc_vars[3];
    frcVID_ = nc_vars[4];
    cellAngleVID_ = nc_vars[5];
    cellLengthVID_ = nc_vars[6];
    timeVID_ = nc_vars[7];
    remd_dimension_ = nc_vars[8];
    indicesVID_ = nc_vars[9];
    ncdebug_ = nc_vars[10];
    ensembleDID_ = nc_vars[11];
    frameDID_ = nc_vars[12];
    atomDID_ = nc_vars[13];
    ncatom_ = nc_vars[14];
    ncatom3_ = nc_vars[15];
    spatialDID_ = nc_vars[16];
    labelDID_ = nc_vars[17];
    cell_spatialDID_ = nc_vars[18];
    cell_angularDID_ = nc_vars[19];
    spatialVID_ = nc_vars[20];
    cell_spatialVID_ = nc_vars[21];
    cell_angularVID_ = nc_vars[22];
  }
}
#endif /* MPI */
#endif /* BINTRAJ */
