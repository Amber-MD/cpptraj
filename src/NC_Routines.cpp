#ifdef BINTRAJ
#include <netcdf.h>
#include "NC_Routines.h"
#include "CpptrajStdio.h"

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
  // Get attr length. Make this fail silently.
  int ncerr = nc_inq_attlen(ncid, vid, attribute, &attlen);
  if (ncerr != NC_NOERR)
    return attrOut;
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
