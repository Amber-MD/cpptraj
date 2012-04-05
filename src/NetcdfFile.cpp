#ifdef BINTRAJ
#  include "netcdf.h"
#endif
#include "NetcdfFile.h"
#include "CpptrajStdio.h"

// DEFINES
#ifdef BINTRAJ
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_LENGTHS "cell_lengths"
#define NCCELL_ANGULAR "cell_angular"
#define NCCELL_ANGLES "cell_angles"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5
#define NCREMD_DIMENSION "remd_dimension"
#define NCREMD_GROUPNUM "remd_groupnum"
#define NCREMD_DIMTYPE "remd_dimtype"
#define NCREMD_INDICES "remd_indices"
#endif

// CONSTRUCTOR
NetcdfFile::NetcdfFile() :
  ncdebug_(0),
  ncid_(-1),
  atomDID_(-1),
  ncatom_(-1),
  ncatom3_(-1),
  coordVID_(-1),
  velocityVID_(-1),
  cellAngleVID_(-1),
  cellLengthVID_(-1),
  spatialDID_(-1),
  labelDID_(-1),
  cell_spatialDID_(-1),
  cell_angularDID_(-1),
  spatialVID_(-1),
  timeVID_(-1),
  cell_spatialVID_(-1),
  cell_angularVID_(-1),
  TempVID_(-1)
{
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 0;
  count_[1] = 0;
  count_[2] = 0;
}

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions() {
  NCTYPE nctype = NC_UNKNOWN;
#ifdef BINTRAJ
  std::string attrText = GetAttrText(NC_GLOBAL, "Conventions");
  if (attrText == "AMBER")
    nctype = NC_AMBERTRAJ;
  else if (attrText == "AMBERRESTART")
    nctype = NC_AMBERRESTART;
  else if (attrText.empty()) 
    mprinterr("Error: Could not get conventions from Netcdf file.\n");
  else {
    mprinterr("Error: Netcdf file: Unrecognized conventions \"%s\".\n",
              attrText.c_str());
    mprinterr("       Expected \"AMBER\" or \"AMBERRESTART\".\n");
  }
#else
   mprintf("Error: Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#endif
   return nctype;
}

#ifdef BINTRAJ
// NetcdfFile::GetAttrText()
/** Get the information about a netcdf attribute with given vid and 
  * attribute text.
  * Since there is no guarantee that NULL char at the end of retrieved string
  * append one.
  */
std::string NetcdfFile::GetAttrText(int vid, const char *attribute) {
  size_t attlen;
  std::string attrOut;
  // Get attr length
  if ( checkNCerr(nc_inq_attlen(ncid_, vid, attribute, &attlen)) ) {
    mprinterr("Error: Getting length for attribute %s\n",attribute); 
    return attrOut;
  }
  // Allocate space for attr text, plus one for NULL char
  char *attrText = new char[ (attlen + 1) ];
  // Get attr text
  if ( checkNCerr(nc_get_att_text(ncid_, vid, attribute, attrText)) ) {
    mprinterr("Error: Getting attribute text for %s\n",attribute);
    delete[] attrText;
    return attrOut;
  }
  // Append NULL char - NECESSARY?
  attrText[attlen]='\0';
  attrOut.assign( attrText );

  return attrOut;
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

// NetcdfFile::SetupCoordinates()
/** Setup ncatom, ncatom3, atomDID, coordVID, spatialDID, spatialVID. 
  * Check units and spatial dimensions.
  */
int NetcdfFile::SetupCoordinates() {
  int spatial;
  atomDID_ = GetDimInfo(NCATOM, &ncatom_);
  if (atomDID_==-1) return 1;
  ncatom3_ = ncatom_ * 3;
  // Get coord info
  if (checkNCerr(nc_inq_varid(ncid_, NCCOORDS, &coordVID_))) {
    mprinterr("Error: Getting coordinate ID\n");
    return 1;
  }
  std::string attrText = GetAttrText(coordVID_, "units");
  if (attrText!="angstrom")
    mprintf("WARNING: Netcdf file has length units of %s - expected angstrom.\n",
            attrText.c_str());
  // Get spatial info
  spatialDID_ = GetDimInfo(NCSPATIAL, &spatial);
  if (spatialDID_==-1) return 1;
  if (spatial!=3) {
    mprinterr("Error: ncOpen: Expected 3 spatial dimenions, got %i\n",spatial);
    return 1;
  }
  if ( checkNCerr(nc_inq_varid(ncid_, NCSPATIAL, &spatialVID_)) ) {
    mprinterr("Error: Getting spatial VID\n");
    return 1;
  }
  return 0;
}

int NetcdfFile::SetupVelocity() {
  if ( nc_inq_varid(ncid_,NCVELO,&velocityVID_)==NC_NOERR ) {
    if (ncdebug_>0) mprintf("    Netcdf file has velocities.\n");
    return 0;
  } 
  velocityVID_=-1;
  return 1;
}

int NetcdfFile::SetupTime() {
  if ( checkNCerr( nc_inq_varid(ncid_, NCTIME, &timeVID_)) ) {
    mprinterr("Getting Netcdf time VID.\n");
    return 1;
  }
  std::string attrText = GetAttrText(timeVID_, "units");
  if (attrText!="picosecond")
    mprintf("Warning: Netcdf file has time units of %s - expected picosecond.\n",
            attrText.c_str());
  return 0;
}

int NetcdfFile::SetupTemperature() {
  if ( nc_inq_varid(ncid_,NCTEMPERATURE,&TempVID_) == NC_NOERR ) {
    if (ncdebug_>0) mprintf("    Netcdf file has replica temperatures.\n");
    return 0;
  } 
  TempVID_=-1;
  return 1;
}

/** \return 0 on success, 1 on error, -1 for no box coords. */
int NetcdfFile::SetupBox(double *boxAngle) {
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
    count_[0]=1; 
    count_[1]=3; 
    count_[2]=0;
    if ( checkNCerr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, boxAngle)) ) {
      mprinterr("Error: Getting cell lengths.\n");
      return 1;
    }
    if ( checkNCerr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, boxAngle)) ) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    return 0;
  }
  // No box information
  return -1;
}

// NetcdfFile::checkNCerr()
bool NetcdfFile::checkNCerr(int ncerr) {
  if ( ncerr != NC_NOERR ) {
    mprinterr("NETCDF Error: %s\n",nc_strerror(ncerr));
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

// NetcdfFile::NC_open()
int NetcdfFile::NC_open(const char* Name) {
  if ( checkNCerr( nc_open( Name, NC_NOWRITE, &ncid_ ) ) )
    return 1;
  return 0;
}

// NetcdfFile::NC_create()
int NetcdfFile::NC_create(const char* Name) {
  if ( checkNCerr( nc_create( Name, NC_64BIT_OFFSET, &ncid_) ) )
    return 1;
  return 0;
}

#endif
