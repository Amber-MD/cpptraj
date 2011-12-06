// NetcdfRoutines
#include "CpptrajStdio.h"
#include <cstdlib>
#ifdef BINTRAJ
#include <cstdio>
#include <cstdarg>
#include "netcdf.h"
#include "NetcdfRoutines.h"

// NetcdfDebug()
/** For use in printing various attributes of a previously opened netcdf file.
  */
void NetcdfDebug(int ncid) {
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  int err,i;
  char *varname;
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
  varname=(char*) malloc(1024*sizeof(char));
  mprintf("========== BEG. NETCDF DEBUG ==========\n");
  err=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  mprintf("nc_inq returned %i\n",err);
  if (err==NC_NOERR)
    mprintf("ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            ndimsp,nvarsp,ngattsp,unlimdimidp);
  else
    mprintf("NETCDF Error occurred.\n");
  // Print name of each variable defined in netcdf file
  mprintf("NC VARIABLES:\n");
  for (i=0; i<nvarsp; i++) {
    err=nc_inq_varname(ncid,i,varname);
    mprintf("  Var %i - ",i);
    if (err==NC_NOERR)
      mprintf("%s\n",varname);
    else
      mprintf("NETCDF Error occured.\n");
  }
  mprintf("==========  END NETCDF DEBUG ==========\n");
  free(varname);
  return;
}

// checkNCerr()
/** Check the given netcdf error status. Print error message and return 1
  * if an error has occurred, otherwise return 0; 
  */
int checkNCerr(int err, const char *message, ...) {
  va_list args;

  if (err!=NC_NOERR) {
    va_start(args,message);
    mprinterr("NETCDF Error (%s): ",nc_strerror(err));
    vfprintf(stderr,message,args);
    mprinterr("\n");
    va_end(args);
    return 1;
  }
  return 0;
}

// GetDimInfo()
/** Return the dimension ID of a given attribute in netcdf file ncid.
  * Also set dimension length.
  */
int GetDimInfo(int ncid, const char *attribute, int *length) {
  int dimID;
  size_t slength;

  *length = 0;
  slength = 0;

  // Get dimid 
  if ( checkNCerr(nc_inq_dimid(ncid, attribute, &dimID),
       "ncGetDimInfo: Getting dimID for attribute %s",attribute)!=0 ) return -1;

  // get Dim length 
  if ( checkNCerr(nc_inq_dimlen(ncid, dimID, &slength),
       "ncGetDimInfo: Getting length for attribute %s",attribute)!=0) return -1;

  *length = (int) slength;
  return dimID;
}

// GetAttrText()
/** Get the information about a netcdf attribute with given vid and 
  * attribute text.
  * Since there is no guarantee that NULL char at the end of retrieved string
  * append one.
  */
char *GetAttrText(int ncid, int vid, const char *attribute) {
  size_t attlen;
  char *attrText;
  // Get attr length
  if ( checkNCerr(nc_inq_attlen(ncid, vid, attribute, &attlen),
       "ncGetAttrText: Getting length for attribute %s\n",attribute)) return NULL;
  // Allocate space for attr text, plus one for NULL char
  attrText = (char*) malloc( (attlen + 1) * sizeof(char));
  // Get attr text
  if ( checkNCerr(nc_get_att_text(ncid, vid, attribute, attrText),
       "ncGetAttrText: Getting attribute text for %s\n",attribute)) {
    free(attrText);
    return NULL;
  }
  // Append NULL char
  attrText[attlen]='\0';

  return attrText;
}
#endif

// GetNetcdfConventions()
/** This is called from CpptrajFile to determine whether filename is a netcdf
  * trajectory or restart file. Return conventions string.
  * Return NULL if error or netcdf not compiled in.
  */
char *GetNetcdfConventions(char *filename) {
#ifdef BINTRAJ
  int ncid;
  char *attrText;
  if (checkNCerr(nc_open(filename,NC_NOWRITE,&ncid),
          "Opening Netcdf file %s for reading",filename)!=0) return NULL;  
  attrText = GetAttrText(ncid,NC_GLOBAL, "Conventions");
  nc_close(ncid);
  return attrText;
#else
   mprintf("Error: Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
   return NULL;
#endif
}

// FloatToDouble()
/** Convert float coords to double coords
  * NOTE: N needs to match up with size of Coord!
  */
void FloatToDouble(double *X, float *Coord, int N) {
  for (int i=0; i<N; i++)
    X[i]=(double) Coord[i];
}

// DoubleToFloat()
/** Convert double coords to float coords
  */
void DoubleToFloat(float *  Coord, double *  X, int N) {
  for (int i=0; i<N; i++)
    Coord[i]=(float) X[i];
}

