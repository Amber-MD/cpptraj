#include <cstdio> // FILE, fopen, fclose
#include "NetcdfFile.h"
# include "CpptrajStdio.h"
#ifdef BINTRAJ
# include <netcdf.h>
# include <cstring> // strlen
# include "NC_Routines.h"
# include "Constants.h"
# include "Version.h"
# include "Frame.h"
# ifdef MPI
#  include "ParallelNetcdf.h"
# endif
# ifdef HAS_HDF5
#  include <limits> // for integer compression
# endif
#endif

// NetcdfFile::GetNetcdfConventions()
/** First check the base format type to determine NetCDF3 vs NetCDF/HDF5.
  * Then check that the file has the proper conventions.
  */
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(NC_FMT_TYPE& btype, const char* fname)
{
  btype = NC_NOTNC;
  NCTYPE nctype = NC_UNKNOWN;
  // Determine base type via magic number
  FILE* infile = fopen(fname, "rb");
  if (infile == 0) return nctype;
  unsigned char buf[8];
  unsigned int nread = fread(buf, sizeof(char), 8, infile);
  //mprintf("DEBUG: '%s' nread=%u %x %x %x %x %x %x %x %x\n", fname, nread, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
  fclose(infile);
  if (nread > 3 && buf[0] == 'C' && buf[1] == 'D' && buf[2] == 'F') {
#   ifdef BINTRAJ
    btype = NC_V3;
#   else
    mprintf("Warning: File '%s' appears to be NetCDF but cpptraj was compiled without NetCDF support.\n", fname);
    return nctype;
#   endif
  } else if (nread > 7 && buf[0] == 0x89 && buf[1] == 0x48 && buf[2] == 0x44 && buf[3] == 0x46 &&
                          buf[4] == 0x0d && buf[5] == 0x0a && buf[6] == 0x1a && buf[7] == 0x0a)
  {
#   ifdef HAS_HDF5
    btype = NC_V4;
#   else
    mprintf("Warning: File '%s' appears to be NetCDF4/HDF5 but cpptraj was compiled without HDF5 support.\n", fname);
    return nctype;
#   endif
  }
# ifdef BINTRAJ
  // NOTE: Do not use checkNCerr so this fails silently. Allows routine to
  //       be used in file autodetection.
  int myNcid;
  if ( nc_open( fname, NC_NOWRITE, &myNcid ) != NC_NOERR )
    return NC_UNKNOWN;
  nctype = GetNetcdfConventions(myNcid);
  nc_close( myNcid ); 
# else
  mprinterr("Error: Compiled without NetCDF support. Recompile with -DBINTRAJ\n");
# endif
  return nctype;
}

/** \return netcdf format type if the file magic number is recognized.
  */

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
#define NCREMDVALUES "remd_values"
#define NCCOMPPOS "compressedpos"
#define NCCOMPPOW "compressedpow"

// CONSTRUCTOR
NetcdfFile::NetcdfFile() :
# ifdef HAS_HDF5
  deflateLevels_((unsigned int)NVID, 0),
  compressedPosVID_(-1),
  compressedFac_(0),
  fchunkSize_(1),
# endif
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
  ensembleSize_(0),
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
  cell_angularVID_(-1),
  RemdValuesVID_(-1)
{
  start_[0] = 0;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 0;
  count_[1] = 0;
  count_[2] = 0;
}

const char* NetcdfFile::ConventionsStr_[] = {
  "AMBER",         // NC_AMBERTRAJ
  "AMBERRESTART",  // NC_AMBERRESTART
  "AMBERENSEMBLE", // NC_AMBERENSEMBLE
  0                // UNKNOWN
};

// NetcdfFile::GetNetcdfConventions()
NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(int ncidIn) {
  NCTYPE nctype = NC_UNKNOWN;
  std::string attrText = NC::GetAttrText(ncidIn, "Conventions");
  if (attrText.empty()) {
    mprinterr("Error: Could not get conventions from NetCDF file.\n");
  } else {
    for (int i = 0; i < (int)NC_UNKNOWN; i++) {
      if (attrText.compare( ConventionsStr_[i] ) == 0) {
        nctype = (NCTYPE)i;
        break;
      }
    }
    if (nctype == NC_UNKNOWN) {
      mprinterr("Error: NetCDF file has unrecognized conventions \"%s\".\n",
                attrText.c_str());
      mprinterr("Error: Expected one of");
      for (int i = 0; i < (int)NC_UNKNOWN; i++)
        mprintf(" \"%s\"", ConventionsStr_[i]);
      mprinterr("\n");
    }
  }
  return nctype;
}

// NetcdfFile::CheckConventionsVersion()
void NetcdfFile::CheckConventionsVersion() {
  std::string attrText = NC::GetAttrText(ncid_, "ConventionVersion");
  if ( attrText != "1.0")
    mprintf("Warning: NetCDF file has ConventionVersion that is not 1.0 (%s)\n", attrText.c_str());
}

/** \return true if temperature VID is defined or a replica dimension is temperature.
  */
bool NetcdfFile::HasTemperatures() const {
  if (TempVID_ != -1)
    return true;
  else if (!remValType_.empty()) {
    for (int idx = 0; idx != remValType_.Ndims(); idx++)
      if (remValType_.DimType(idx) == ReplicaDimArray::TEMPERATURE)
        return true;
  }
  return false;
}

/** \return true if a replica dimension is pH. TODO put inside ReplicaDimArray
  */
bool NetcdfFile::Has_pH() const {
  if (!remValType_.empty()) {
    for (int idx = 0; idx != remValType_.Ndims(); idx++)
      if (remValType_.DimType(idx) == ReplicaDimArray::PH)
        return true;
  }
  return false;
}

/** \return true if a replica dimension is RedOx.
  */
bool NetcdfFile::HasRedOx() const {
  if (!remValType_.empty()) {
    for (int idx = 0; idx != remValType_.Ndims(); idx++)
      if (remValType_.DimType(idx) == ReplicaDimArray::REDOX)
        return true;
  }
  return false;
}

// NetcdfFile::SetupFrameDim()
/** Get the frame dimension ID and # of frames (ncframe). */
int NetcdfFile::SetupFrameDim() {
  frameDID_ = NC::GetDimInfo( ncid_, NCFRAME, ncframe_ );
  if (frameDID_==-1) return 1;
  return 0;
}

/** Get the ensemble dimension ID and size. */
int NetcdfFile::SetupEnsembleDim() {
  ensembleSize_ = 0;
  ensembleDID_ = NC::GetDimInfo( ncid_, NCENSEMBLE, ensembleSize_ );
  if (ensembleDID_ == -1) return 0;
  return ensembleSize_;
}

// NetcdfFile::SetupCoordsVelo()
/** Setup ncatom, ncatom3, atomDID, coordVID, spatialDID, spatialVID,
  * velocityVID, frcVID. Check units and spatial dimensions.
  */
int NetcdfFile::SetupCoordsVelo(bool useVelAsCoords, bool useFrcAsCoords) {
  if (useVelAsCoords && useFrcAsCoords) {
    mprinterr("Error: Cannot use both velocities and forces as coords - specify one only.\n");
    return 1;
  }
  // Get atoms info
  atomDID_ = NC::GetDimInfo( ncid_, NCATOM, ncatom_ );
  if (atomDID_==-1) return 1;
  ncatom3_ = ncatom_ * 3;
  // Check for integer compression - requires HDF5
  int localCompressedPosVID = -1;
  if ( nc_inq_varid(ncid_, NCCOMPPOS, &localCompressedPosVID) == NC_NOERR ) {
    if (ncdebug_ > 0) mprintf("\tNetCDF file has integer-compressed coordinates.\n");
  }
  if (localCompressedPosVID != -1) {
#   ifdef HAS_HDF5
    int compressedPowVID = -1;
    compressedPosVID_ = localCompressedPosVID;
    // Get compressed power var ID
    if ( nc_inq_varid(ncid_, NCCOMPPOW, &compressedPowVID) != NC_NOERR ) {
      mprinterr("Error: NetCDF file has integer-compressed coordinates but no compress power variable ID.\n");
      return 1;
    }
    // Get compressed power
    int power = 0;
    if (NC::CheckErr(nc_get_var_int(ncid_, compressedPowVID, &power))) {
      mprinterr("Error: Reading compressed power factor.\n");
      return 1;
    }
    // Calculate compressed factor
    if (calcCompressFactor(power)) return 1;
    // Allocate temporary space for integer array
    itmp_.resize( Ncatom3() );
#   else
    mprinterr("Error: Integer-compressed NetCDF trajectories requires cpptraj compiled with HDF5 support.\n");
    return 1;
#   endif /* HAS_HDF5 */
  }
  // Get coord info
  coordVID_ = -1;
  if ( nc_inq_varid(ncid_, NCCOORDS, &coordVID_) == NC_NOERR ) {
    if (ncdebug_ > 0) mprintf("\tNetCDF file has coordinates.\n");
    std::string attrText = NC::GetAttrText(ncid_, coordVID_, "units");
    if (attrText!="angstrom")
      mprintf("Warning: NetCDF file has length units of %s - expected angstrom.\n",
              attrText.c_str());
  }
  // Get spatial info
  int spatial;
  spatialDID_ = NC::GetDimInfo( ncid_, NCSPATIAL, spatial );
  if (spatialDID_==-1) return 1;
  if (spatial!=3) {
    mprinterr("Error: Expected 3 spatial dimensions, got %i\n",spatial);
    return 1;
  }
  if ( NC::CheckErr(nc_inq_varid(ncid_, NCSPATIAL, &spatialVID_)) ) {
    mprintf("Warning: Could not get spatial VID. File may not be Amber NetCDF compliant.\n");
    mprintf("Warning: Assuming spatial variables are 'x', 'y', 'z'\n");
  } else {
    start_[0] = 0;
    count_[0] = 3;
    char xyz[3];
    if (NC::CheckErr(nc_get_vara_text(ncid_, spatialVID_, start_, count_, xyz))) {
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
    if (ncdebug_>0) mprintf("\tNetCDF file has velocities.\n");
  }
  // Get force info
  frcVID_ = -1;
  if ( nc_inq_varid(ncid_, NCFRC, &frcVID_) == NC_NOERR ) {
    if (ncdebug_>0) mprintf("\tNetCDF file has forces.\n");
  }
  // Return a error if no coords, velocities, or forces
  if ( coordVID_ == -1 && velocityVID_ == -1 && frcVID_ == -1 ) {
    mprinterr("Error: NetCDF file has no coordinates, velocities, or forces.\n");
    return 1;
  }
  // If using velocities/forces as coordinates, swap them now.
  if (useVelAsCoords) {
    if (velocityVID_ == -1) {
      mprinterr("Error: Cannot use velocities as coordinates; no velocities present.\n");
      return 1;
    }
    mprintf("\tUsing velocities as coordinates.\n");
    coordVID_ = velocityVID_;
    velocityVID_ = -1;
  } else if (useFrcAsCoords) {
    if (frcVID_ == -1) {
      mprinterr("Error: Cannot use forces as coordinates; no forces present.\n");
      return 1;
    }
    mprintf("\tUsing forces as coordinates.\n");
    coordVID_ = frcVID_;
    frcVID_ = -1;
  }
  // Get overall replica and coordinate indices
  crdidxVID_ = -1;
  if ( nc_inq_varid(ncid_, NCREMD_REPIDX, &repidxVID_) == NC_NOERR ) {
      if (ncdebug_>0) mprintf("\tNetCDF file has overall replica indices.\n");
    if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_CRDIDX, &crdidxVID_)) ) {
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
    std::string attrText = NC::GetAttrText(ncid_, timeVID_, "units");
    if (attrText!="picosecond")
      mprintf("Warning: NetCDF file has time units of %s - expected picosecond.\n",
              attrText.c_str());
    // Check for time values which have NOT been filled, which was possible
    // with netcdf trajectories created by older versions of ptraj/cpptraj.
    if (ncframe_ > 0 && myType_ == NC_AMBERTRAJ) {
      float time;
      start_[0] = 0; count_[0] = 1;
      if (NC::CheckErr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time))) {
        mprinterr("Error: Getting time value for NetCDF file.\n");
        return -1;
      }
      if (time == NC_FILL_FLOAT) {
        mprintf("Warning: NetCDF file time variable defined but empty. Disabling.\n");
        timeVID_ = -1;
      } else {
        // If first 2 values are 0, this is another indication of a bad time variable.
        if (ncframe_ > 1) {
          float time1;
          start_[0] = 1;
          if (NC::CheckErr(nc_get_vara_float(ncid_, timeVID_, start_, count_, &time1))) {
            mprinterr("Error: Getting second time value for NetCDF file.\n");
            return -1;
          }
          if (time1 == 0 && time1 == time)
          {
            mprintf("Warning: NetCDF file time variable defined but all zero. Disabling.\n");
            timeVID_ = -1;
          }
        }
      }
    }
    return 0;
  }
  timeVID_=-1;
  return 1;
}

// NetcdfFile::SetupTemperature()
/** Determine if Netcdf file contains temperature; set up temperature VID. */
void NetcdfFile::SetupTemperature() {
  TempVID_ = -1;
  if ( nc_inq_varid(ncid_, NCTEMPERATURE, &TempVID_) == NC_NOERR ) {
    if (ncdebug_ > 0) mprintf("\tNetCDF file has replica temperatures.\n");
  }
}

// NetcdfFile::SetupMultiD()
/** Determine if Netcdf file contains multi-D REMD info. If so set the
  * number of replica dimensions (remd_dimension_) and figure out
  * the dimension types (remdDim)
  */
int NetcdfFile::SetupMultiD() {
  int dimensionDID;
  remd_dimension_ = 0;
  if ( nc_inq_dimid(ncid_, NCREMD_DIMENSION, &dimensionDID) == NC_NOERR)
  {
    // Although this is a second call to dimid, makes for easier code
    if ( (dimensionDID = NC::GetDimInfo(ncid_, NCREMD_DIMENSION, remd_dimension_))==-1 )
      return -1;
    if (ncdebug_ > 0)
      mprintf("\tNetCDF file has multi-D REMD info, %i dimensions.\n",remd_dimension_);
    // Ensure valid # dimensions
    if (remd_dimension_ < 1) {
      if (ncdebug_ > 0)
        mprintf("\tNumber of REMD dimensions is less than 1.\n");
      remd_dimension_ = 0;
    } else {
      // Start and count for groupnum and dimtype, allocate mem
      start_[0]=0; 
      start_[1]=0; 
      start_[2]=0;
      count_[0]=remd_dimension_; 
      count_[1]=0; 
      count_[2]=0;
      std::vector<int> remd_dimtype( remd_dimension_ );
      // Get dimension types
      int dimtypeVID;
      if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_DIMTYPE, &dimtypeVID)) ) {
        mprinterr("Error: Getting dimension type variable ID for each dimension.\n");
        return -1;
      }
      if ( NC::CheckErr(nc_get_vara_int(ncid_, dimtypeVID, start_, count_, &remd_dimtype[0])) ) {
        mprinterr("Error: Getting dimension type in each dimension.\n");
        return -1;
      }
      // Get VID for replica indices
      if ( NC::CheckErr(nc_inq_varid(ncid_, NCREMD_INDICES, &indicesVID_)) ) {
        mprinterr("Error: Getting replica indices variable ID.\n");
        return -1;
      }
      // Store type of each dimension TODO should just have one netcdf coordinfo
      remDimType_.clear();
      for (int dim = 0; dim < remd_dimension_; ++dim)
        remDimType_.AddRemdDimension( remd_dimtype[dim] );
    }
  }

  // Get VID for replica values
  if ( nc_inq_varid(ncid_, NCREMDVALUES, &RemdValuesVID_) == NC_NOERR ) {
    if (ncdebug_ > 0) mprintf("\tNetCDF file has replica values.\n");
    remValType_.clear();
    if (remd_dimension_ > 0) {
      remValType_ = remDimType_;
    } else {
      // Probably 1D
      int ndims = 0;
      if (NC::CheckErr(nc_inq_varndims(ncid_, RemdValuesVID_, &ndims))) {
        mprinterr("Error: Checking number of dimensions for REMD_VALUES.\n");
        return 1;
      }
      if (ndims < 2) {
        mprintf("Warning: New style (Amber >= V18) remd values detected with < 2 dimensions.\n"
                "Warning:  Assuming temperature.\n");
        remValType_.AddRemdDimension( ReplicaDimArray::TEMPERATURE );
      } else {
        mprinterr("Error: New style (Amber >= V18) remd values detected with > 1 dimension\n"
                  "Error:   but no multi-D replica info present.\n");
        return 1;
      }
    }
    RemdValues_.assign( remValType_.Ndims(), 0 );
  }

  return 0; 
}

// NetcdfFile::SetupBox()
/** \return 0 on success, 1 on error, -1 for no box coords. */
int NetcdfFile::SetupBox() {
  nc_box_.SetNoBox();
  if ( nc_inq_varid(ncid_, NCCELL_LENGTHS, &cellLengthVID_) == NC_NOERR ) {
    if (NC::CheckErr( nc_inq_varid(ncid_, NCCELL_ANGLES, &cellAngleVID_) )) {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (ncdebug_ > 0) mprintf("\tNetCDF Box information found.\n");
    // If present, get box angles. These will be used to determine the box 
    // type in TrajectoryFile.
    start_[0]=0; 
    start_[1]=0; 
    start_[2]=0;
    start_[3]=0;
    switch (myType_) {
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
    double boxCrd[6]; /// XYZ ABG
    if ( NC::CheckErr(nc_get_vara_double(ncid_, cellLengthVID_, start_, count_, boxCrd )) )
    {
      mprinterr("Error: Getting cell lengths.\n");
      return 1;
    }
    if ( NC::CheckErr(nc_get_vara_double(ncid_, cellAngleVID_, start_, count_, boxCrd+3)) )
    {
      mprinterr("Error: Getting cell angles.\n");
      return 1;
    }
    if (ncdebug_ > 0) mprintf("\tNetCDF Box: XYZ={%f %f %f} ABG={%f %f %f}\n",
                              boxCrd[0], boxCrd[1], boxCrd[2], boxCrd[3], boxCrd[4], boxCrd[5]);
    if (nc_box_.SetupFromXyzAbg( boxCrd )) {
      mprintf("Warning: NetCDF file unit cell variables appear to be empty; disabling box.\n");
      cellLengthVID_ = -1;
      cellAngleVID_ = -1;
      nc_box_.SetNoBox();
      return -1;
    }
    return 0;
  }
  // No box information
  return -1;
}

/** Read REMD related values. */
int NetcdfFile::ReadRemdValues(Frame& frm) {
  // FIXME assuming start_ is set
  count_[0] = 1; // 1 frame
  if ( repidxVID_ != -1)
    nc_get_vara_int(ncid_, repidxVID_, start_, count_, frm.repidxPtr());
  if ( crdidxVID_ != -1)
    nc_get_vara_int(ncid_, crdidxVID_, start_, count_, frm.crdidxPtr());
  if ( RemdValuesVID_ != -1 ) {
    count_[1] = remd_dimension_; // # dimensions
    if ( NC::CheckErr(nc_get_vara_double(ncid_, RemdValuesVID_, start_, count_, &RemdValues_[0])) )
    {
      mprinterr("Error: Getting replica values\n");
      return 1;
    }
    for (int idx = 0; idx != remValType_.Ndims(); ++idx)
    {
      if (remValType_.DimType(idx) == ReplicaDimArray::TEMPERATURE) {
        frm.SetTemperature( RemdValues_[idx] );
        //mprintf("DEBUG: T= %g\n", frm.Temperature());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::PH) {
        frm.Set_pH( RemdValues_[idx] );
        //mprintf("DEBUG: pH= %g\n", frm.pH());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::REDOX) {
        frm.SetRedOx( RemdValues_[idx] );
        //mprintf("DEBUG: RedOx= %g\n", frm.RedOx());
      }
    }
  }
  return 0;
}

/** Set up a NetCDF file for reading. */
int NetcdfFile::NC_setupRead(std::string const& fname, NCTYPE expectedType, int expectedNatoms,
                             bool useVelAsCoords, bool useFrcAsCoords, int debugIn)
{
  ncdebug_ = debugIn;
  // If file is open, close it.
  if (ncid_ != -1) NC_close();
  // Open read
  if (NC_openRead( fname.c_str() ) != 0) {
    mprinterr("Error: Could not open NetCDF file '%s' for read setup.\n", fname.c_str());
    return 1;
  }
  // Sanity check
  myType_ = GetNetcdfConventions(ncid_);
  if ( myType_ != expectedType ) {
    mprinterr("Error: NetCDF file conventions do not include \"%s\"\n",
              ConventionsStr_[expectedType]);
    return 1;
  }
  // This will warn if conventions are not 1.0 
  CheckConventionsVersion();
  // Get the title
  nctitle_ = NC::GetAttrText(ncid_, "title");
  // Get frame info if necessary.
  if (myType_ == NC_AMBERTRAJ || myType_ == NC_AMBERENSEMBLE) {
    if (SetupFrameDim() != 0) return 1;
    if (Ncframe() < 1) {
      mprinterr("Error: NetCDF file has no frames.\n");
      return 1;
    }
    // Get ensemble info if necessary
    if (myType_ == NC_AMBERENSEMBLE) {
      if (SetupEnsembleDim() < 1) {
        mprinterr("Error: Could not get ensemble dimension info.\n");
        return 1;
      }
    }
  }
  // Setup atom-dimension-related variables. 
  if ( SetupCoordsVelo( useVelAsCoords, useFrcAsCoords ) != 0 ) return 1;
  // Check that specified number of atoms matches expected number.
  if (Ncatom() != expectedNatoms) {
    mprinterr("Error: Number of atoms in NetCDF file (%i) does not match number\n"
              "Error:  in associated topology (%i)!\n", Ncatom(), expectedNatoms);
    return 1; 
  }
  // Setup Time - FIXME: Allowed to fail silently
  SetupTime();
  // Box info
  if (SetupBox() == 1) // 1 indicates an error
    return 1;
  // Replica Temperatures - FIXME: Allowed to fail silently
  SetupTemperature();
  // Replica Dimensions
  if ( SetupMultiD() == -1 ) return 1;
  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;
  if (ncdebug_ > 1) NC::Debug(ncid_);
  NC_close();
  return 0;
}
  

/** \return Coordinate info corresponding to current setup. */
CoordinateInfo NetcdfFile::NC_coordInfo() const {
  // TODO the 'false' is for step info. Enable this in the future when time
  //      is present.
  return CoordinateInfo( ensembleSize_, remDimType_, nc_box_,
                         HasCoords(), HasVelocities(), HasForces(), 
                         HasTemperatures(), Has_pH(), HasRedOx(),
                         HasTimes(), false, (repidxVID_ != -1), (crdidxVID_ != -1),
                         (RemdValuesVID_ != -1) );
}

// NetcdfFile::NC_openRead()
int NetcdfFile::NC_openRead(std::string const& Name) {
  if (Name.empty()) return 1;
  if ( NC::CheckErr( nc_open( Name.c_str(), NC_NOWRITE, &ncid_ ) ) )
    return 1;
  return 0;
}

// NetcdfFile::NC_close()
void NetcdfFile::NC_close() {
  if (ncid_ == -1) return;
  bool err = NC::CheckErr( nc_close(ncid_) );
  if (ncdebug_ > 0 && !err)
    mprintf("Successfully closed ncid %i\n",ncid_);
  ncid_ = -1;
}

// =============================================================================
// NetcdfFile::NC_openWrite()
int NetcdfFile::NC_openWrite(std::string const& Name) {
  if (Name.empty()) return 1;
  if ( NC::CheckErr( nc_open( Name.c_str(), NC_WRITE, &ncid_ ) ) )
    return 1;
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
  if ( NC::CheckErr( nc_redef( ncid_ ) ) ) return 1;
  // Define eptot, bins, temp0
  if ( NC::CheckErr( nc_def_var(ncid_, NCEPTOT, NC_DOUBLE, 1, dimensionID, &eptotVID)) ) {
    mprinterr("Error: defining eptot variable ID.\n");
    return 1;
  }
  if (hasBins) {
    if ( NC::CheckErr( nc_def_var(ncid_, NCBINS, NC_INT, 1, dimensionID, &binsVID)) ) {
      mprinterr("Error: defining bins variable ID.\n");
      return 1;
    }
  } else
    binsVID = -1;
  if (NC_defineTemperature(dimensionID, 0)) return 1;
  // Random seed, make global
  if ( NC::CheckErr( nc_put_att_int(ncid_, NC_GLOBAL, "iseed", NC_INT, 1, &iseed) ) ) {
    mprinterr("Error: setting random seed attribute.\n");
    return 1;
  }
  // End definitions
  if (NC::CheckErr(nc_enddef(ncid_))) {
    mprinterr("NetCDF error on ending definitions.");
    return 1;
  }
  // Write temperature
  if (NC::CheckErr(nc_put_var_double(ncid_,TempVID_,&reservoirT)) ) {
    mprinterr("Error: Writing reservoir temperature.\n");
    return 1;
  }
  return 0;
}

/** Set remDimDID appropriate for given type. */
void NetcdfFile::SetRemDimDID(int remDimDID, int* dimensionID) const {
  switch (myType_) {
    case NC_AMBERENSEMBLE:
      dimensionID[0] = frameDID_;
      dimensionID[1] = ensembleDID_;
      dimensionID[2] = remDimDID;
      break;
    case NC_AMBERTRAJ:
      dimensionID[0] = frameDID_;
      dimensionID[1] = remDimDID;
      break;
    case NC_AMBERRESTART:
      dimensionID[0] = remDimDID;
      break;
    default:
      mprinterr("Internal Error: SetRemDimDID(): Unrecognized type.\n");
  }
}

/** Variable ID descriptions. */
const char* NetcdfFile::vidDesc_[NVID] = {
  "coordinates",     // V_COORDS
  "velocities",      // V_VEL
  "forces",          // V_FRC
  "temp0",           // V_TEMP
  "cell_lengths",    // V_BOXL
  "cell_angles",     // V_BOXA
  "time",            // V_TIME
  "remd_indices",    // V_IND
  "remd_repidx",     // V_RIDX
  "remd_crdidx",     // V_CIDX
  "remd_values"      // V_REMDVALS
};

/** Dimension ID descriptions. */
const char* NetcdfFile::didDesc_[NDID] = {
  "frame",           // D_FRAME
  "atom",            // D_ATOM
  "spatial"          // D_SPATIAL
};

// -----------------------------------------------------------------------------
// This section contains HDF5-related functionality

/** Set variable compression level if supported. */
int NetcdfFile::NC_setDeflate(VidType vtype, int varid) const
{
# ifdef HAS_HDF5
  if (deflateLevels_[vtype] > 0) {
    // TODO shuffle integer types?
    if ( NC::CheckErr( nc_def_var_deflate(ncid_, varid, 0, 1, deflateLevels_[vtype]) ) ) {
      mprinterr("Error: Setting compression for '%s' variable.\n", vidDesc_[vtype]);
      return 1;
    }
  }
# else
  mprintf("Warning: Setting NetCDF variable compression requires compiling with HDF5 support.\n");
# endif
  return 0;
}

/** Increase frame chunk size for variable if supported. */
int NetcdfFile::NC_setFrameChunkSize(VidType vtype, int varid) const
{
# ifdef HAS_HDF5
  if (fchunkSize_ > 1) {
    // Get number of dimensions
    int ndims = 0;
    if ( NC::CheckErr( nc_inq_varndims(ncid_, varid, &ndims) ) ) {
      mprinterr("Error: getting # dims for '%s' variable.\n", vidDesc_[vtype]);
      return 1;
    }
    if (ndims < 1) {
      mprintf("Warning: NC_setFrameChunkSize: Variable '%s' has no dimensions.\n",
              vidDesc_[vtype]);
      return 0;
    }
    // Get dimension IDs
    std::vector<int> dimids(ndims);
    if ( NC::CheckErr( nc_inq_var(ncid_, varid, 0, 0, 0, &dimids[0], 0) ) ) {
      mprinterr("Error: getting dimension IDs for '%s' variable.\n", vidDesc_[vtype]);
      return 1;
    }
    // Allocate space for chunk sizes 
    std::vector<size_t> chunkSizes( ndims );
    // Set frame chunk size
    int err = NC_setVarDimChunkSizes(vtype, varid, fchunkSize_, dimids, frameDID_, chunkSizes);
    return err;
  }
# else
  mprintf("Warning: Setting NetCDF frame chunk size requires compiling with HDF5 support.\n");
# endif
  return 0;
}

# ifdef HAS_HDF5
/** Set desired compression level for specified variable if supported. */
int NetcdfFile::SetCompression(VidType vtype, int deflateLevelIn) {
  if (ncdebug_ > 0)
    mprintf("DEBUG: Setting compression for VIDTYPE %i to %i\n", (int)vtype, deflateLevelIn);
  deflateLevels_[vtype] = deflateLevelIn;
  return 0;
}

/** Set desired compression level for all variables if supported. */
int NetcdfFile::SetCompression(int deflateLevelIn) {
  mprintf("\tSetting NetCDF variable compression level to %i\n", deflateLevelIn);
  int err = 0;
  for (int i = 0; i != (int)NVID; i++)
    err += SetCompression( (VidType)i, deflateLevelIn );
  return err;
}

/* Set desired frame chunk size if supported. */
int NetcdfFile::SetFrameChunkSize(int fchunkSizeIn) {
  mprintf("\tSetting frame chunk size to %i\n", fchunkSizeIn);
  fchunkSize_ = fchunkSizeIn;
  return 0;
}

/** If target dimID is -1 multiply existing chunk sizes for variable
  * by chunkFac. Otherwise multiple chunk size matching target dimID
  * only.
  */
int NetcdfFile::NC_setVarDimChunkSizes(VidType vtype, int varid, int chunkFac,
                                       std::vector<int> const& dimids, int tgtDimId,
                                       std::vector<size_t>& chunkSizes)
const
{
  // Sanity checks.
  if (dimids.empty() || chunkSizes.empty()) {
    mprinterr("Internal Error: NC_setVarDimChunkSizes called with empty arrays.\n");
    return 1;
  }
  // Get chunk sizes and storage type
  int storagep = 0;
  if ( NC::CheckErr( nc_inq_var_chunking(ncid_, varid, &storagep, &chunkSizes[0]) ) ) {
    mprinterr("Error: getting existing chunk sizes for '%s' variable.\n", vidDesc_[vtype]);
    return 1;
  }
  if (storagep != NC_CHUNKED) {
    mprinterr("Internal Error: NC_setVarDimChunkSizes called for VID that is not chunked '%s'\n",
              vidDesc_[vtype]);
    return 1;
  }
  // Loop over dimension chunk sizes
  for (unsigned int dim = 0; dim != dimids.size(); dim++) {
    if (tgtDimId == -1 || tgtDimId == dimids[dim]) {
      mprintf("DEBUG: '%s' dim %u chunk size = %zu new size= %i\n",
              vidDesc_[vtype], dim, chunkSizes[dim], chunkSizes[dim]*chunkFac);
      // Set new chunk size
      chunkSizes[dim] *= chunkFac;
    }
  }
  if ( NC::CheckErr( nc_def_var_chunking(ncid_, varid, NC_CHUNKED, &chunkSizes[0]) ) ) {
    mprinterr("Error: Setting chunk sizes for '%s' variable.\n", vidDesc_[vtype]);
    return 1;
  }
  return 0;
}
 
/** Set compressedFac_ to given power of 10 (min 1). */
int NetcdfFile::calcCompressFactor(int power) {
  if (power < 1) {
    mprinterr("Internal Error: calcCompressFactor called with power of 10 < 1\n");
    return 1;
  }
  compressedFac_ = 10.0;
  for (int i = 1; i < power; i++)
    compressedFac_ *= 10.0;
  mprintf("\tConverting floats to integer using factor: x%g\n", compressedFac_);
  return 0;
}

/** Prepare trajectory for adding coords converted to integer. */
int NetcdfFile::NC_createIntCompressed(int power)
{
  if (ncid_ == -1) {
    mprinterr("Internal Error: NC_createCompressed() called with ncid -1\n");
    return 1;
  }
  if (frameDID_ == -1 || atomDID_ == -1 || spatialDID_ == -1) {
    mprinterr("Internal Error: NC_createCompressed() called before NC_create().\n");
    return 1;
  }
  // Warn about low precision powers
  if (power < 1) {
    mprinterr("Error: Integer compression power < 1 not allowed.\n");
    return 1;
  } else if (power < 4) {
    mprintf("Warning: Using extremely low precision for integer compression.\n"
            "Warning: Energy error will be on the order of 2E-%i kcal/mol/atom\n", power);
    mprintf("Warning: Consider using integer power >= 4.\n");
  } else
    mprintf("Warning: Using lossy compression.\n"
            "Warning: Energy error will be on the order of 2E-%i kcal/mol/atom\n", power);
  int dimensionID[NC_MAX_VAR_DIMS];
  // Place file back in define mode
  if (NC::CheckErr( nc_redef( ncid_ ) )) return 1;
  // Define variable to hold converted coords
  dimensionID[0] = frameDID_;
  dimensionID[1] = atomDID_;
  dimensionID[2] = spatialDID_;
  if (NC::CheckErr( nc_def_var(ncid_, NCCOMPPOS, NC_INT, 3, dimensionID, &compressedPosVID_) )) {
    mprinterr("Error: defining compressed positions VID.\n");
    return 1;
  }
  // Set compression for converted coords. Use 1 if not set.
  if (deflateLevels_[V_COORDS] == 0) {
    mprintf("Warning: Using default compression level for coords.\n");
    deflateLevels_[V_COORDS] = 1;
  }
  // TODO should this be an option?
  int ishuffle = 1;
  if (ncdebug_ > 0) {
    if (ishuffle == 0)
      mprintf("DEBUG: Integer shuffle is off.\n");
    else
      mprintf("DEBUG: Integer shuffle is on.\n");
  }
  if ( NC::CheckErr( nc_def_var_deflate(ncid_, compressedPosVID_, ishuffle, 1, deflateLevels_[V_COORDS]) ) ) {
    mprinterr("Error: Setting compression level %i for integer compressed coords.\n", deflateLevels_[V_COORDS]);
    return 1;
  }
  if (NC_setFrameChunkSize(V_COORDS, compressedPosVID_)) return 1;
  // Define variable to hold conversion power
  int compressedPowVID = -1;
  if (NC::CheckErr( nc_def_var(ncid_, NCCOMPPOW, NC_INT, 0, dimensionID, &compressedPowVID) )) {
    mprinterr("Error: defining compressed power factor VID.\n");
    return 1;
  }
  // End definitions
  if (NC::CheckErr(nc_enddef(ncid_))) {
    mprinterr("NetCDF error on ending definitions (NC_createCompressed).");
    return 1;
  }
  // Calculate compressed factor; write compressed power of 10
  if (calcCompressFactor(power)) return 1;
  if (NC::CheckErr(nc_put_var_int(ncid_, compressedPowVID, &power)) ) {
    mprinterr("Error: Writing compressed power factor.\n");
    return 1;
  }
  // Allocate temporary space for integer array
  itmp_.resize( Ncatom3() );

  return 0;
}

/** Write integer-compressed coords. */
int NetcdfFile::NC_writeIntCompressed(Frame const& frmOut) {
  // Convert to integer
  long int maxval = (long int)std::numeric_limits<int>::max();
  for (int idx = 0; idx != frmOut.size(); idx++)
  {
    // Try some overflow protection
    long int ii = (long int)(frmOut[idx] * compressedFac_);
    if (ii > maxval || ii < -maxval) {
      mprinterr("Error: Coordinate %i frame %i (%g) is too large to convert to int.\n",
                idx+1, ncframe_+1, frmOut[idx]);
      mprinterr("Error: A smaller integer compression factor must be used.\n");
      return 1;
    }
    itmp_[idx] = (int)ii;
    //itmp_[idx] = (int)(frmOut[idx] * compressedFac_);
  }
  //mprintf("DEBUG: atom 0 xyz={ %20.10f %20.10f %20.10f } ixyz= { %20i %20i %20i }\n",
  //        frmOut[0], frmOut[1], frmOut[2], itmp_[0], itmp_[1], itmp_[2]);
  //  Write array
  start_[0] = ncframe_;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;

  if (NC::CheckErr(nc_put_vara_int(ncid_, compressedPosVID_, start_, count_, &itmp_[0]))) {
    mprinterr("Error: NetCDF writing compressed coordinates frame %i\n", ncframe_+1);
    return 1;
  }
  return 0;
}

/** Read integer-compressed coords. */
int NetcdfFile::NC_readIntCompressed(int set, Frame& frmIn) {
  // Read array
  start_[0] = set;
  start_[1] = 0;
  start_[2] = 0;
  count_[0] = 1;
  count_[1] = Ncatom();
  count_[2] = 3;
  if (NC::CheckErr(nc_get_vara_int(ncid_, compressedPosVID_, start_, count_, &itmp_[0]))) {
    mprinterr("Error: NetCDF reading compressed coordinates frame %i\n", set+1);
    return 1;
  }
  // Convert from integer
  for (int idx = 0; idx != frmIn.size(); idx++)
    frmIn[idx] = (double)(itmp_[idx]) / compressedFac_; // TODO convert to 1/fac first?
  return 0;
}
// -----------------------------------------------------------------------------
#endif /* HAS_HDF5 */

// NetcdfFile::NC_defineTemperature()
int NetcdfFile::NC_defineTemperature(int* dimensionID, int NDIM) {
  if (NC::CheckErr(nc_def_var(ncid_,NCTEMPERATURE,NC_DOUBLE,NDIM,dimensionID,&TempVID_))) {
    mprinterr("NetCDF error on defining temperature.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,TempVID_,"units",6,"kelvin"))) {
    mprinterr("NetCDF error on defining temperature units.\n");
    return 1;
  }
  if (NC_setDeflate(V_TEMP, TempVID_)) return 1;
  return 0;
}

/** Create default NetCDF version 3 file. */
int NetcdfFile::NC_create(std::string const& Name, NCTYPE typeIn, int natomIn,
                          CoordinateInfo const& coordInfo, std::string const& title, int debugIn) 
{
  return (NC_create(NC_V3, Name, typeIn, natomIn, coordInfo, title, debugIn));
}


// NetcdfFile::NC_create()
int NetcdfFile::NC_create(NC_FMT_TYPE wtypeIn, std::string const& Name, NCTYPE typeIn, int natomIn,
                          CoordinateInfo const& coordInfo, std::string const& title, int debugIn) 
{
  if (Name.empty()) return 1;
  int dimensionID[NC_MAX_VAR_DIMS];
  int NDIM;
  nc_type dataType;
  myType_ = typeIn;
  ncdebug_ = debugIn;

  if (ncdebug_>1)
    mprintf("DEBUG: NC_create: '%s'  natom=%i  %s\n",
            Name.c_str(),natomIn, coordInfo.InfoString().c_str());
  int cmode;
  if (wtypeIn == NC_V3)
    cmode = NC_64BIT_OFFSET;
  else if (wtypeIn == NC_V4)
    cmode = NC_NETCDF4;
  else {
    mprinterr("Internal Error: Unspecified base format type given in NC_create().\n");
    return 1;
  }
  if ( NC::CheckErr( nc_create( Name.c_str(), cmode, &ncid_) ) )
    return 1;

  ncatom_ = natomIn;
  ncatom3_ = ncatom_ * 3;
  
  // Set number of dimensions based on file type
  switch (myType_) {
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
      mprinterr("Error: NC_create (%s): Unrecognized type (%i)\n",Name.c_str(),(int)myType_);
      return 1;
  }

  // Ensemble dimension for ensemble
  if (myType_ == NC_AMBERENSEMBLE) {
    if (coordInfo.EnsembleSize() < 1) {
      mprinterr("Internal Error: NetcdfFile: ensembleSize < 1\n");
      return 1;
    }
    if ( NC::CheckErr(nc_def_dim(ncid_, NCENSEMBLE, coordInfo.EnsembleSize(), &ensembleDID_)) ) {
      mprinterr("Error: Defining ensemble dimension.\n");
      return 1;
    }
    dimensionID[1] = ensembleDID_;
  }
  // Frame dimension for traj
  ncframe_ = 0;
  if (myType_ == NC_AMBERTRAJ || myType_ == NC_AMBERENSEMBLE) {
    if ( NC::CheckErr( nc_def_dim( ncid_, NCFRAME, NC_UNLIMITED, &frameDID_)) ) {
      mprinterr("Error: Defining frame dimension.\n");
      return 1;
    }
    // Since frame is UNLIMITED, it must be lowest dim.
    dimensionID[0] = frameDID_;
  }
  // Time variable and units
  if (coordInfo.HasTime()) {
    if ( NC::CheckErr( nc_def_var(ncid_, NCTIME, dataType, NDIM-2, dimensionID, &timeVID_)) ) {
      mprinterr("Error: Defining time variable.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text(ncid_, timeVID_, "units", 10, "picosecond")) ) {
      mprinterr("Error: Writing time VID units.\n");
      return 1;
    }
    if (NC_setDeflate(V_TIME, timeVID_)) return 1;
  }
  // Spatial dimension and variable
  if ( NC::CheckErr( nc_def_dim( ncid_, NCSPATIAL, 3, &spatialDID_)) ) {
    mprinterr("Error: Defining spatial dimension.\n");
    return 1;
  }
  dimensionID[0] = spatialDID_;
  if ( NC::CheckErr( nc_def_var( ncid_, NCSPATIAL, NC_CHAR, 1, dimensionID, &spatialVID_)) ) {
    mprinterr("Error: Defining spatial variable.\n"); 
    return 1;
  }
  // Atom dimension
  if ( NC::CheckErr( nc_def_dim( ncid_, NCATOM, ncatom_, &atomDID_)) ) {
    mprinterr("Error: Defining atom dimension.\n");
    return 1;
  }
  // Setup dimensions for Coords/Velocity
  // NOTE: THIS MUST BE MODIFIED IF NEW TYPES ADDED
  switch (myType_) {
    case NC_AMBERENSEMBLE:
      dimensionID[0] = frameDID_;
      dimensionID[1] = ensembleDID_;
      dimensionID[2] = atomDID_;
      dimensionID[3] = spatialDID_;
      break;
    case NC_AMBERTRAJ:
      dimensionID[0] = frameDID_;
      dimensionID[1] = atomDID_;
      dimensionID[2] = spatialDID_;
      break;
    case NC_AMBERRESTART:
      dimensionID[0] = atomDID_;
      dimensionID[1] = spatialDID_;
      break;
    case NC_UNKNOWN:
      mprinterr("Internal Error: Unknown type passed to NC_create()\n");
      return 1;
  }
  // Coord variable
  if (coordInfo.HasCrd()) {
    if ( NC::CheckErr( nc_def_var( ncid_, NCCOORDS, dataType, NDIM, dimensionID, &coordVID_)) ) {
      mprinterr("Error: Defining coordinates variable.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text( ncid_, coordVID_, "units", 8, "angstrom")) ) {
      mprinterr("Error: Writing coordinates variable units.\n");
      return 1;
    }
    if (NC_setDeflate(V_COORDS, coordVID_)) return 1;
    if (NC_setFrameChunkSize(V_COORDS, coordVID_)) return 1;
  }
  // Velocity variable
  if (coordInfo.HasVel()) {
    if ( NC::CheckErr( nc_def_var( ncid_, NCVELO, dataType, NDIM, dimensionID, &velocityVID_)) ) {
      mprinterr("Error: Defining velocities variable.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text( ncid_, velocityVID_, "units", 19, "angstrom/picosecond")) )
    {
      mprinterr("Error: Writing velocities variable units.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_double( ncid_, velocityVID_, "scale_factor", NC_DOUBLE, 1, 
                                        &Constants::AMBERTIME_TO_PS)) )
    {
      mprinterr("Error: Writing velocities scale factor.\n");
      return 1;
    }
    if (NC_setDeflate(V_VEL, velocityVID_)) return 1;
    if (NC_setFrameChunkSize(V_VEL, velocityVID_)) return 1;
  }
  // Force variable
  if (coordInfo.HasForce()) {
    if ( NC::CheckErr( nc_def_var( ncid_, NCFRC, dataType, NDIM, dimensionID, &frcVID_)) ) {
      mprinterr("Error: Defining forces variable\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text( ncid_, frcVID_, "units", 25, "kilocalorie/mole/angstrom")) )
    {
      mprinterr("Error: Writing forces variable units.\n");
      return 1;
    }
    if (NC_setDeflate(V_FRC, frcVID_)) return 1;
    if (NC_setFrameChunkSize(V_FRC, frcVID_)) return 1;
  }
  // Replica Temperature
  if (coordInfo.HasTemp() && !coordInfo.UseRemdValues()) {
    // NOTE: Setting dimensionID should be OK for Restart, will not be used.
    dimensionID[0] = frameDID_;
    if ( NC_defineTemperature( dimensionID, NDIM-2 ) ) return 1;
  }
  // Overall replica index
  if (coordInfo.HasRepIdx()) {
    dimensionID[0] = frameDID_;
    if (NC::CheckErr(nc_def_var(ncid_, NCREMD_REPIDX, NC_INT, NDIM-2, dimensionID, &repidxVID_)))
    {
      mprinterr("Error: Defining replica idx variable ID.\n");
      return 1;
    }
    if (NC_setDeflate(V_RIDX, repidxVID_)) return 1;
    //if (NC_setFrameChunkSize(V_RIDX, repidxVID_)) return 1;
  }
  // Overall coordinate index
  if (coordInfo.HasCrdIdx()) {
    dimensionID[0] = frameDID_;
    if (NC::CheckErr(nc_def_var(ncid_, NCREMD_CRDIDX, NC_INT, NDIM-2, dimensionID, &crdidxVID_)))
    {
      mprinterr("Error: Defining coordinate idx variable ID.\n");
      return 1;
    }
    if (NC_setDeflate(V_CIDX, crdidxVID_)) return 1;
    //if (NC_setFrameChunkSize(V_CIDX, crdidxVID_)) return 1;
  }
  // Replica indices
  int remDimTypeVID = -1;
  int remDimDID = -1;
  if (coordInfo.HasReplicaDims()) {
    // Define number of replica dimensions
    remd_dimension_ = coordInfo.ReplicaDimensions().Ndims();
    if ( NC::CheckErr(nc_def_dim(ncid_, NCREMD_DIMENSION, remd_dimension_, &remDimDID)) ) {
      mprinterr("Error: Defining replica indices dimension.\n");
      return 1;
    }
    dimensionID[0] = remDimDID;
    // For each dimension, store the type
    if ( NC::CheckErr(nc_def_var(ncid_, NCREMD_DIMTYPE, NC_INT, 1, dimensionID, &remDimTypeVID)) ) 
    {
      mprinterr("Error: Defining replica dimension type variable.\n");
      return 1;
    }
    // Need to store the indices of replica in each dimension each frame
    SetRemDimDID(remDimDID, dimensionID);
    if (NC::CheckErr(nc_def_var(ncid_, NCREMD_INDICES, NC_INT, NDIM-1, dimensionID, &indicesVID_)))
    {
      mprinterr("Error: Defining replica indices variable ID.\n");
      return 1;
    }
    if (NC_setDeflate(V_IND, indicesVID_)) return 1;
    if (NC_setFrameChunkSize(V_IND, indicesVID_)) return 1;
    // TODO: Determine if groups are really necessary for restarts. If not, 
    // remove from AmberNetcdf.F90.
  }
  // Replica values
  if (coordInfo.UseRemdValues()) {
    remValType_.clear();
    if (coordInfo.HasReplicaDims()) {
      SetRemDimDID(remDimDID, dimensionID);
      if (NC::CheckErr(nc_def_var(ncid_, NCREMDVALUES, NC_DOUBLE, NDIM-1,
                                  dimensionID, &RemdValuesVID_)))
      {
        mprinterr("Error: defining replica values variable ID.\n");
        return 1;
      }
      RemdValues_.resize( remd_dimension_ );
      for (int idx = 0; idx != coordInfo.ReplicaDimensions().Ndims(); idx++)
        remValType_.AddRemdDimension( coordInfo.ReplicaDimensions().DimType(idx) );
    } else {
      dimensionID[0] = frameDID_;
      if (NC::CheckErr(nc_def_var(ncid_, NCREMDVALUES, NC_DOUBLE, NDIM-2,
                                  dimensionID, &RemdValuesVID_)))
      {
        mprinterr("Error: defining replica values variable ID.\n");
        return 1;
      }
      RemdValues_.resize( 1 );
      // FIXME assuming temperature
      remValType_.AddRemdDimension( ReplicaDimArray::TEMPERATURE );
    }
    if (NC_setDeflate(V_REMDVALS, RemdValuesVID_)) return 1;
    //if (NC_setFrameChunkSize(V_REMDVALS, RemdValuesVID_)) return 1;
  }
  // Box Info
  if (coordInfo.HasBox()) {
    // Check x-aligned
    if (!coordInfo.TrajBox().Is_X_Aligned())
      mprintf("Warning: Unit cell is not X-aligned. Box cannot be properly stored as Amber NetCDF\n");
    // Cell Spatial
    if ( NC::CheckErr( nc_def_dim( ncid_, NCCELL_SPATIAL, 3, &cell_spatialDID_)) ) {
      mprinterr("Error: Defining cell spatial dimension.\n");
      return 1;
    }
    dimensionID[0] = cell_spatialDID_;
    if ( NC::CheckErr( nc_def_var(ncid_, NCCELL_SPATIAL, NC_CHAR, 1, dimensionID, &cell_spatialVID_)))
    {
      mprinterr("Error: Defining cell spatial variable.\n");
      return 1;
    }
    // Cell angular
    if ( NC::CheckErr( nc_def_dim( ncid_, NCLABEL, NCLABELLEN, &labelDID_)) ) {
      mprinterr("Error: Defining label dimension.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_def_dim( ncid_, NCCELL_ANGULAR, 3, &cell_angularDID_)) ) {
      mprinterr("Error: Defining cell angular dimension.\n"); 
      return 1;
    }
    dimensionID[0] = cell_angularDID_;
    dimensionID[1] = labelDID_;
    if ( NC::CheckErr( nc_def_var( ncid_, NCCELL_ANGULAR, NC_CHAR, 2, dimensionID, 
                                 &cell_angularVID_)) )
    {
      mprinterr("Error: Defining cell angular variable.\n");
      return 1;
    }
    // Setup dimensions for Box
    // NOTE: This must be modified if more types added
    int boxdim;
    if (myType_ == NC_AMBERENSEMBLE) {
      dimensionID[0] = frameDID_;
      dimensionID[1] = ensembleDID_;
      boxdim = 2;
    } else if (myType_ == NC_AMBERTRAJ) {
      dimensionID[0] = frameDID_;
      boxdim = 1;
    } else {
      boxdim = 0;
    }
    dimensionID[boxdim] = cell_spatialDID_;
    if ( NC::CheckErr( nc_def_var( ncid_, NCCELL_LENGTHS, NC_DOUBLE, NDIM-1, dimensionID,
                                 &cellLengthVID_)) )
    {
      mprinterr("Error: Defining cell length variable.\n"); 
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text( ncid_, cellLengthVID_, "units", 8, "angstrom")) ) {
      mprinterr("Error: Writing cell length variable units.\n");
      return 1;
    }
    dimensionID[boxdim] = cell_angularDID_;
    if ( NC::CheckErr( nc_def_var( ncid_, NCCELL_ANGLES, NC_DOUBLE, NDIM-1, dimensionID,
                                 &cellAngleVID_)) )
    {
      mprinterr("Error: Defining cell angle variable.\n");
      return 1;
    }
    if ( NC::CheckErr( nc_put_att_text( ncid_, cellAngleVID_, "units", 6, "degree")) ) {
      mprinterr("Error: Writing cell angle variable units.\n");
      return 1;
    }
    if (NC_setDeflate(V_BOXL, cellLengthVID_)) return 1;
    if (NC_setFrameChunkSize(V_BOXL, cellLengthVID_)) return 1;
    if (NC_setDeflate(V_BOXA, cellAngleVID_)) return 1;
    if (NC_setFrameChunkSize(V_BOXA, cellAngleVID_)) return 1;
  }

  // Attributes
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"title",title.size(),title.c_str())) ) {
    mprinterr("Error: Writing title.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"application",5,"AMBER")) ) {
    mprinterr("Error: Writing application.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"program",7,"cpptraj")) ) {
    mprinterr("Error: Writing program.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"programVersion",
                                   strlen(CPPTRAJ_INTERNAL_VERSION), CPPTRAJ_INTERNAL_VERSION)))
  {
    mprinterr("Error: Writing program version.\n");
    return 1;
  }
  // Write conventions based on type 
  std::string cStr( ConventionsStr_[myType_] );
  if (NC::CheckErr(nc_put_att_text(ncid_, NC_GLOBAL, "Conventions", cStr.size(), cStr.c_str())))
  {
    mprinterr("Error: Writing conventions.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"ConventionVersion",3,"1.0")) ) {
    mprinterr("Error: Writing conventions version.\n");
    return 1;
  }
  
  // Set fill mode
  if (NC::CheckErr(nc_set_fill(ncid_, NC_NOFILL, dimensionID))) {
    mprinterr("Error: NetCDF setting fill value.\n");
    return 1;
  }

  // End netcdf definitions
  if (NC::CheckErr(nc_enddef(ncid_))) {
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
  if (NC::CheckErr(nc_put_vara_text(ncid_, spatialVID_, start_, count_, xyz))) {
    mprinterr("Error on NetCDF output of spatial VID 'x', 'y' and 'z'");
    return 1;
  }
  if ( coordInfo.HasBox() ) {
    xyz[0] = 'a'; 
    xyz[1] = 'b'; 
    xyz[2] = 'c';
    if (NC::CheckErr(nc_put_vara_text(ncid_, cell_spatialVID_, start_, count_, xyz))) {
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
    if (NC::CheckErr(nc_put_vara_text(ncid_, cell_angularVID_, start_, count_, abc))) {
      mprinterr("Error on NetCDF output of cell angular VID 'alpha', 'beta ' and 'gamma'");
      return 1;
    }
  }

  // Store the type of each replica dimension.
  if (coordInfo.HasReplicaDims()) {
    remDimType_.clear();
    ReplicaDimArray const& remdDim = coordInfo.ReplicaDimensions();
    start_[0] = 0;
    count_[0] = remd_dimension_;
    std::vector<int> tempDims( remd_dimension_ );
    for (int i = 0; i < remd_dimension_; ++i) {
      tempDims[i] = remdDim[i];
      remDimType_.AddRemdDimension( remdDim[i] );
    }
    if (NC::CheckErr(nc_put_vara_int(ncid_, remDimTypeVID, start_, count_, &tempDims[0]))) {
      mprinterr("Error: writing replica dimension types.\n");
      return 1;
    }
  }
  if (ncdebug_ > 1) NC::Debug(ncid_);

  return 0;
}

/** Write REMD-related values. */
int NetcdfFile::WriteRemdValues(Frame const& frm) {
  // FIXME assuming start_ is set
  count_[0] = 1; // 1 frame
  if ( repidxVID_ != -1)
    nc_put_vara_int(ncid_, repidxVID_, start_, count_, frm.repidxPtr());
  if ( crdidxVID_ != -1)
    nc_put_vara_int(ncid_, crdidxVID_, start_, count_, frm.crdidxPtr());
  if ( RemdValuesVID_ != -1 ) {
    for (int idx = 0; idx != remValType_.Ndims(); ++idx)
    {
      if (remValType_.DimType(idx) == ReplicaDimArray::TEMPERATURE) {
        RemdValues_[idx] = frm.Temperature();
        //mprintf("DEBUG: T= %g\n", frm.Temperature());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::PH) {
        RemdValues_[idx] = frm.pH();
        //mprintf("DEBUG: pH= %g\n", frm.pH());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::REDOX) {
        RemdValues_[idx] = frm.RedOx();
        //mprintf("DEBUG: RedOx= %g\n", frm.RedOx());
      }
    }
    count_[1] = remd_dimension_; // # dimensions
    if ( NC::CheckErr(nc_put_vara_double(ncid_, RemdValuesVID_, start_, count_, &RemdValues_[0])) )
    {
      mprinterr("Error: Writing replica values\n");
      return 1;
    }
  }
  return 0;
}

#ifdef MPI
#ifdef HAS_PNETCDF
int NetcdfFile::parallelWriteRemdValues(int set, Frame const& frm) {
  MPI_Offset pstart_[2];
  MPI_Offset pcount_[2];
  pstart_[0] = set;
  pstart_[1] = 0;
  pcount_[0] = 1;               // 1 frame
  pcount_[1] = remd_dimension_; // # dimensions 
  // REMD values
  if ( repidxVID_ != -1)
    ncmpi_put_vara_int(ncid_, repidxVID_, pstart_, pcount_, frm.repidxPtr());
  if ( crdidxVID_ != -1)
    ncmpi_put_vara_int(ncid_, crdidxVID_, pstart_, pcount_, frm.crdidxPtr());
  if ( RemdValuesVID_ != -1 ) {
    for (int idx = 0; idx != remValType_.Ndims(); ++idx)
    {
      if (remValType_.DimType(idx) == ReplicaDimArray::TEMPERATURE) {
        RemdValues_[idx] = frm.Temperature();
        //mprintf("DEBUG: T= %g\n", frm.Temperature());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::PH) {
        RemdValues_[idx] = frm.pH();
        //mprintf("DEBUG: pH= %g\n", frm.pH());
      } else if (remValType_.DimType(idx) == ReplicaDimArray::REDOX) {
        RemdValues_[idx] = frm.RedOx();
        //mprintf("DEBUG: RedOx= %g\n", frm.RedOx());
      }
    }
    pcount_[1] = remd_dimension_; // # dimensions
    if ( NC::CheckErr(ncmpi_put_vara_double(ncid_, RemdValuesVID_, pstart_, pcount_, &RemdValues_[0])) )
    {
      mprinterr("Error: Writing replica values\n");
      return 1;
    }
  }
  return 0;
}
#endif
#endif

// =============================================================================
// NetcdfFile::DebugIndices()
void NetcdfFile::DebugIndices() const {
  mprintf("DBG: Start={%zu, %zu, %zu, %zu} Count={%zu, %zu, %zu, %zu}\n",
         start_[0], start_[1], start_[2], start_[3],
         count_[0], count_[1], count_[2], count_[3]);
}

// NetcdfFile::DebugVIDs()
void NetcdfFile::DebugVIDs() const {
  rprintf("TempVID_=%i  coordVID_=%i  velocityVID_=%i frcVID_=%i  cellAngleVID_=%i"
          "  cellLengthVID_=%i  indicesVID_=%i\n",
          TempVID_, coordVID_, velocityVID_, frcVID_, cellAngleVID_, cellLengthVID_, indicesVID_);
}

#ifdef MPI
void NetcdfFile::Sync(Parallel::Comm const& commIn) {
  static const unsigned int NCVARS_SIZE = 29;
  int nc_vars[NCVARS_SIZE];
  if (commIn.Master()) {
    nc_vars[0]  = ncframe_;
    nc_vars[1]  = TempVID_;
    nc_vars[2]  = coordVID_;
    nc_vars[3]  = velocityVID_;
    nc_vars[4]  = frcVID_;
    nc_vars[5]  = cellAngleVID_;
    nc_vars[6]  = cellLengthVID_;
    nc_vars[7]  = timeVID_;
    nc_vars[8]  = remd_dimension_;
    nc_vars[9]  = indicesVID_;
    nc_vars[10] = repidxVID_;
    nc_vars[11] = crdidxVID_;
    nc_vars[12] = ensembleSize_;
    nc_vars[13] = ncdebug_;
    nc_vars[14] = ensembleDID_;
    nc_vars[15] = frameDID_;
    nc_vars[16] = atomDID_;
    nc_vars[17] = ncatom_;
    nc_vars[18] = ncatom3_;
    nc_vars[19] = spatialDID_;
    nc_vars[20] = labelDID_;
    nc_vars[21] = cell_spatialDID_;
    nc_vars[22] = cell_angularDID_;
    nc_vars[23] = spatialVID_;
    nc_vars[24] = cell_spatialVID_;
    nc_vars[25] = cell_angularVID_;
    nc_vars[26] = RemdValuesVID_;
    nc_vars[27] = remDimType_.Ndims();
    nc_vars[28] = remValType_.Ndims();
    commIn.MasterBcast( nc_vars, NCVARS_SIZE, MPI_INT );
  } else {
    // Non-master
    commIn.MasterBcast( nc_vars, NCVARS_SIZE, MPI_INT );
    ncframe_         = nc_vars[0];
    TempVID_         = nc_vars[1];
    coordVID_        = nc_vars[2];
    velocityVID_     = nc_vars[3];
    frcVID_          = nc_vars[4];
    cellAngleVID_    = nc_vars[5];
    cellLengthVID_   = nc_vars[6];
    timeVID_         = nc_vars[7];
    remd_dimension_  = nc_vars[8];
    indicesVID_      = nc_vars[9];
    repidxVID_       = nc_vars[10];
    crdidxVID_       = nc_vars[11];
    ensembleSize_    = nc_vars[12];
    ncdebug_         = nc_vars[13];
    ensembleDID_     = nc_vars[14];
    frameDID_        = nc_vars[15];
    atomDID_         = nc_vars[16];
    ncatom_          = nc_vars[17];
    ncatom3_         = nc_vars[18];
    spatialDID_      = nc_vars[19];
    labelDID_        = nc_vars[20];
    cell_spatialDID_ = nc_vars[21];
    cell_angularDID_ = nc_vars[22];
    spatialVID_      = nc_vars[23];
    cell_spatialVID_ = nc_vars[24];
    cell_angularVID_ = nc_vars[25];
    RemdValuesVID_   = nc_vars[26];
    remDimType_.assign( nc_vars[27], ReplicaDimArray::UNKNOWN );
    remValType_.assign( nc_vars[28], ReplicaDimArray::UNKNOWN );
    RemdValues_.resize( nc_vars[28] );
  }
  if (!remDimType_.empty())
    commIn.MasterBcast( remDimType_.Ptr(), remDimType_.Ndims(), MPI_INT );
  if (!remValType_.empty())
    commIn.MasterBcast( remValType_.Ptr(), remValType_.Ndims(), MPI_INT );
}
#endif /* MPI */
#endif /* BINTRAJ */
