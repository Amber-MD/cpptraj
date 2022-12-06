#include "NC_Routines.h"
#include "CpptrajStdio.h"
#include <cstdio>
#ifdef BINTRAJ
# include <netcdf.h>
#endif

// NC::GetConventions()
NC::ConventionsType NC::GetConventions(std::string const& fname) {
  ConventionsType nctype = NC_UNKNOWN;
# ifdef BINTRAJ
  // NOTE: Do not use checkNCerr so this fails silently. Allows routine to
  //       be used in file autodetection.
  int myNcid;
  if ( nc_open( fname.c_str(), NC_NOWRITE, &myNcid ) != NC_NOERR )
    return NC_UNKNOWN;
  nctype = GetConventions(myNcid);
  nc_close( myNcid ); 
# else
  mprinterr("Error: Compiled without NetCDF support. Recompile with -DBINTRAJ\n");
# endif
  return nctype;
}

// NOTE: Must be kept in sync with FormatType in NC_Routines.h
static const char* NC_FmtTypeStr_[] = {
  "NetCDF3",
  "NetCDF4/HDF5",
  "Not NetCDF"
};

// NC::fmtTypeStr()
const char* NC::fmtTypeStr(FormatType ftype) {
  return NC_FmtTypeStr_[ftype];
}

/** \return Format type of specified file. */
NC::FormatType NC::GetFormatType(std::string const& fnameIn) {
  // Determine base type via magic number
  FILE* infile = fopen(fnameIn.c_str(), "rb");
  if (infile == 0) return NC_NOTNC;
  unsigned char buf[8];
  unsigned int nread = fread(buf, sizeof(char), 8, infile);
  //mprintf("DEBUG: '%s' nread=%u %x %x %x %x %x %x %x %x\n", fname, nread, buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7]);
  fclose(infile);
  if (nread > 3 && buf[0] == 'C' && buf[1] == 'D' && buf[2] == 'F') {
#   ifdef BINTRAJ
    return NC_V3;
#   else
    mprintf("Warning: File '%s' appears to be NetCDF3 but cpptraj was compiled without NetCDF support.\n", fnameIn.c_str());
    return NC_NOTNC;
#   endif
  } else if (nread > 7 && buf[0] == 0x89 && buf[1] == 0x48 && buf[2] == 0x44 && buf[3] == 0x46 &&
                          buf[4] == 0x0d && buf[5] == 0x0a && buf[6] == 0x1a && buf[7] == 0x0a)
  {
#   ifdef HAS_HDF5
    return NC_V4;
#   else
    mprintf("Warning: File '%s' appears to be NetCDF4/HDF5 but cpptraj was compiled without HDF5 support.\n", fnameIn.c_str());
    return NC_NOTNC;
#   endif
  }
  //mprintf("DEBUG: Type: %s\n", NcFmtTypeStr_[btype]);
  return NC_NOTNC;
}

#ifdef BINTRAJ

// NC::CheckErr()
bool NC::CheckErr(int ncerr) {
  if ( ncerr != NC_NOERR ) {
    mprintf("%s\n", nc_strerror(ncerr));
    return true;
  }
  return false;
}

/** Get the information about a NetCDF text attribute with given vid and
  * attribute text. Since there is no guarantee that null char at the
  * end of retrieved string, append one.
  */
std::string NC::GetAttrText(int ncid, int vid, const char* attribute) {
  size_t attlen;
  std::string attrOut;
  // Get attr length
  int ncerr = nc_inq_attlen(ncid, vid, attribute, &attlen);
  if (ncerr != NC_NOERR) {
    // Only print error message if something other than attribute not present happened.
    if (ncerr != NC_ENOTATT)
      CheckErr( ncerr );
    return attrOut;
  }
  //if ( CheckErr(nc_inq_attlen(ncid, vid, attribute, &attlen)) ) {
  //  mprintf("Warning: Getting length for attribute '%s'\n",attribute);
  //  return attrOut;
  //}
  // Check attribute type
  int xtypep = -1;
  if ( CheckErr(nc_inq_atttype(ncid, vid, attribute, &xtypep)) ) {
    mprintf("Warning: Problem getting attribute type for '%s'\n", attribute);
    return attrOut;
  }
  //mprintf("DEBUG: Attribute %s type %i\n", attribute, xtypep);

  // Get attr text
  if ( xtypep == NC_CHAR ) {
    // Allocate space for attr text, plus one for null char
    char *attrText = new char[ (attlen + 1) ];
    if ( CheckErr(nc_get_att_text(ncid, vid, attribute, attrText)) ) {
      mprintf("Warning: Problem getting attribute text for '%s'\n", attribute);
      delete[] attrText;
      return attrOut;
    }
    // Append null char - NECESSARY?
    attrText[attlen]='\0';
    attrOut.assign( attrText );
    delete[] attrText;
  } else if ( xtypep == NC_STRING ) {
#   ifdef HAS_HDF5
    if (attlen > 1) {
      mprinterr("Error: Variable attribute %s has %zu strings, expected only 1.\n",
                attribute, attlen);
      return attrOut;
    }
    char* attrString;
    if ( CheckErr(nc_get_att_string(ncid, vid, attribute, &attrString)) ) {
      mprintf("Warning: Problem getting attribute string for '%s'\n", attribute);
      return attrOut;
    }
    attrOut.assign( attrString );
    nc_free_string(attlen, &attrString);
#   else
    mprinterr("Internal Error: Attribute type 'string' only supported with NetCDF4/HDF5.\n"
              "Internal Error: Recompile with Netcdf4/HDF5 support.\n");
#   endif
  } else {
    mprinterr("Error: Attribute %s has unhandled type.\n", attribute);
    return attrOut;
  }

  return attrOut;
}

// NC::GetAttrText()
std::string NC::GetAttrText(int ncid, const char* attribute) {
  return GetAttrText(ncid, NC_GLOBAL, attribute);
}

// NC::GetDimInfo()
int NC::GetDimInfo(int ncid, const char* attribute, unsigned int& length) {
  int dimID;
  size_t slength = 0;

  length = 0;
  // Get dimid 
  if ( CheckErr(nc_inq_dimid(ncid, attribute, &dimID)) ) {
    mprinterr("Error: Getting dimID for attribute %s\n",attribute);
    return -1;
  }
  // get Dim length 
  if ( CheckErr(nc_inq_dimlen(ncid, dimID, &slength)) ) {
    mprinterr("Error: Getting length for attribute %s\n",attribute);
    return -1;
  }
  length = (unsigned int)slength;
  return dimID;
}

// FIXME For backwards compat. only
int NC::GetDimInfo(int ncid, const char* attribute, int& length) {
  unsigned int ulen;
  int dimID = GetDimInfo(ncid, attribute, ulen);
  length = (int)ulen;
  return dimID;
}

// NC::Debug()
void NC::Debug(int ncid) {
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  char varname[NC_MAX_NAME+1];
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
  int err = nc_inq(ncid, &ndimsp, &nvarsp, &ngattsp, &unlimdimidp);
  mprintf("nc_inq returned %i\n", err);
  if (err == NC_NOERR)
    mprintf("ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            ndimsp, nvarsp, ngattsp, unlimdimidp);
  else
    mprintf("NETCDF Error occurred.\n");
  // Print name of each dimension defined in netcdf file
  mprintf("NC DIMS:\n");
  for (int i = 0; i < ndimsp; i++) {
    err = nc_inq_dimname(ncid, i, varname);
    mprintf("  Dim %i - ", i);
    if (err == NC_NOERR)
      mprintf("%s\n", varname);
    else
      mprintf("NETCDF Error occured.\n");
  }
  // Print name of each variable defined in netcdf file
  mprintf("NC VARIABLES:\n");
  for (int i = 0; i < nvarsp; i++) {
    err = nc_inq_varname(ncid, i, varname);
    mprintf("  Var %i - ", i);
    if (err == NC_NOERR)
      mprintf("%s\n", varname);
    else
      mprintf("NETCDF Error occured.\n");
  }
  mprintf("==========  END NETCDF DEBUG ==========\n");
}

/// NOTE: Must be kept in sync with NC::ConventionsType
static const char* NC_ConventionsStr_[] = {
  "AMBER",           // NC_AMBERTRAJ
  "AMBERRESTART",    // NC_AMBERRESTART
  "AMBERENSEMBLE",   // NC_AMBERENSEMBLE
  "CPPTRAJ_CMATRIX", // NC_CPPTRAJCMATRIX
  "CPPTRAJ_DATA",    // NC_CPPTRAJDATA
  0                  // UNKNOWN
};

// NC::conventionsStr()
const char* NC::conventionsStr(ConventionsType nctype) {
  return NC_ConventionsStr_[nctype];
}

// NC::GetConventions()
/** Get conventions. No warning/error messages if no conventions or no
  * recognized conventions found.
  */
NC::ConventionsType NC::GetConventions(int ncid) {
  return GetConventions(ncid, false);
}

// NC::GetConventions()
NC::ConventionsType NC::GetConventions(int ncidIn, bool verbose) {
  ConventionsType nctype = NC_UNKNOWN;
  std::string attrText = GetAttrText(ncidIn, "Conventions");
  if (attrText.empty()) {
    if (verbose)
      mprintf("Warning: Could not get conventions from NetCDF file.\n");
  } else {
    for (int i = 0; i < (int)NC_UNKNOWN; i++) {
      if (attrText.compare( NC_ConventionsStr_[i] ) == 0) {
        nctype = (ConventionsType)i;
        break;
      }
    }
    if (verbose && nctype == NC_UNKNOWN) {
      mprintf("Warning: NetCDF file has unrecognized conventions \"%s\".\n",
              attrText.c_str());
      mprintf("Warning: Expected one of");
      for (int i = 0; i < (int)NC_UNKNOWN; i++)
        mprintf(" \"%s\"", NC_ConventionsStr_[i]);
      mprintf("\n");
    }
  }
  return nctype;
}

/** Add Conventions string to given NetCDF id. */
int NC::PutConventions(int ncid, ConventionsType nctype) {
  std::string cStr( conventionsStr(nctype) );
  if (CheckErr(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", cStr.size(), cStr.c_str())))
  {
    mprinterr("Error: Writing NetCDF Conventions for '%s'.\n", cStr.c_str());
    return 1;
  }
  return 0;
}

// ----- NetCDF4/HDF5 routines -------------------------------------------------
#ifdef HAS_HDF5
/** \return Array containing group names. Also set array with corresponding ncids. */
std::vector<std::string> NC::GetGroupNames(int ncid, std::vector<int>& NcidArray) {
  std::vector<std::string> GroupNames;
  int numgrps;
  nc_inq_grps( ncid, &numgrps, NULL );
  if (numgrps < 1)
    return GroupNames;

  //mprintf("DEBUG: Netcdf file contains %i groups.\n", numgrps);
  GroupNames.reserve( numgrps );
  NcidArray.assign( numgrps, -1 );
  int* ncids = &NcidArray[0];
  nc_inq_grps( ncid, NULL, ncids );
  for (int ii = 0; ii < numgrps; ii++) {
    //mprintf("\tncid %i", ncids[ii]);
    size_t gnamelen;
    nc_inq_grpname_len( ncids[ii], &gnamelen );
    char* gname = new char[ gnamelen+1 ];
    nc_inq_grpname( ncids[ii], gname );
    //mprintf(" %s\n", gname);
    GroupNames.push_back( gname );
    delete[] gname;
  }
  return GroupNames;
}

/** \return Array containing group names. */
std::vector<std::string> NC::GetGroupNames(int ncid) {
  std::vector<int> NcidArray;
  return GetGroupNames(ncid, NcidArray);
}
#endif /* HAS_HDF5 */
// -----------------------------------------------------------------------------
#endif /* BINTRAJ */
