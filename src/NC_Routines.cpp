#ifdef BINTRAJ
#include <netcdf.h>
#include "NC_Routines.h"
#include "CpptrajStdio.h"

bool NC::CheckErr(int ncerr) {
  if ( ncerr != NC_NOERR ) {
    mprintf("%s\n", nc_strerror(ncerr));
    return true;
  }
  return false;
}

std::string NC::GetAttrText(int ncid, int vid, const char* attribute) {
  size_t attlen;
  std::string attrOut;
  // Get attr length
  if ( CheckErr(nc_inq_attlen(ncid, vid, attribute, &attlen)) ) {
    mprintf("Warning: Getting length for attribute '%s'\n",attribute);
    return attrOut;
  }
  // Allocate space for attr text, plus one for null char
  char *attrText = new char[ (attlen + 1) ];
  // Get attr text
  if ( CheckErr(nc_get_att_text(ncid, vid, attribute, attrText)) ) {
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

std::string NC::GetAttrText(int ncid, const char* attribute) {
  return GetAttrText(ncid, NC_GLOBAL, attribute);
}

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
#endif
