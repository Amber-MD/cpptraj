#ifdef HASNETCDF
/* 
 * This file contains a collection of routines designed for reading
 * (and writing?) netcdf files used with amber.
 * Dan Roe 10-2008
 * Original implementation of netcdf in Amber by Jon Mongan.
 */
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring> // For title length
#include "netcdf.h"
#include "AmberNetcdf.h"

#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_ANGULAR "cell_angular"
#define NCCOORDS "coordinates"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5

// CONSTRUCTOR
AmberNetcdf::AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Constructor\n");
  ncid=-1;
  frameID=-1;
  ncframe=-1;
  atomID=-1;
  ncatom=-1;
  Coord=NULL;
  coordID=-1;
  cellAngleID=-1;
  cellLengthID=-1;

  spatialID=-1;
  labelID=-1;
  cell_spatialID=-1;
  cell_angularID=-1;
  spatialVID=-1;
  timeVID=-1;
  cell_spatialVID=-1;
  cell_angularVID=-1;
  TempVarID=-1;

  Coord=NULL;
} 

AmberNetcdf::~AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  if (Coord!=NULL) free(Coord);
  // NOTE: Need to close file?
}

/* DAN DEBUG 
 * dan_netcdf_debug()
 * For use in printing various attributes of a previously opened netcdf file.
 */
void AmberNetcdf::Debug() {
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  int err,i;
  char *varname;

  /* ncid:    NetCDF ID, from a previous call to nc open or nc create.
   * ndimsp:  Pointer to location for returned number of dimensions defined for 
   *         this netCDF dataset.
   * nvarsp:  Pointer to location for returned number of variables defined for 
   *         this netCDF dataset.
   * ngattsp: Pointer to location for returned number of global attributes 
   *         defined for this netCDF dataset.
   * unlimdimidp: 
   *  Pointer to location for returned ID of the unlimited dimension, if 
   *  there is one for this netCDF dataset. If no unlimited length 
   *  dimension has been defined, -1 is returned.
   */
  varname=(char*) malloc(1024*sizeof(char));
  fprintf(stdout,"========== BEG. NETCDF DEBUG ==========\n");
  err=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  fprintf(stdout,"nc_inq returned %i\n",err);
  if (err==NC_NOERR)
    fprintf(stdout,"ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
            ndimsp,nvarsp,ngattsp,unlimdimidp);
  else
    fprintf(stdout,"NETCDF Error occurred.\n");
  /* Print name of each variable defined in netcdf file */
  fprintf(stdout,"NC VARIABLES:\n");
  for (i=0; i<nvarsp; i++) {
    err=nc_inq_varname(ncid,i,varname);
    fprintf(stdout,"  Var %i - ",i);
    if (err==NC_NOERR)
      fprintf(stdout,"%s\n",varname);
    else
      fprintf(stdout,"NETCDF Error occured.\n");
  }

  fprintf(stdout,"==========  END NETCDF DEBUG ==========\n");

  free(varname);
  return;
}


/*
 * checkNCerr()
 * NOTE: Should err be a class member, ncerr?
 */
int AmberNetcdf::checkNCerr(int err, const char *message, ...) {
  va_list args;

  if (err!=NC_NOERR) {
    va_start(args,message);
    fprintf(stderr,"NETCDF Error (%s): ",nc_strerror(err));
    vfprintf(stderr,message,args);
    fprintf(stderr,"\n");
    va_end(args);
    return 1;
  }
  return 0;
}
    

/*
 * ncGetDimInfo()
 * Return the dimension ID of a given attribute in netcdf file ncid.
 * Also set dimension length.
 */ 
int AmberNetcdf::GetDimInfo(const char *attribute, int *length) {
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

/*
 * AmberNetcdf::close()
 */
void AmberNetcdf::close() {
  checkNCerr(nc_close(ncid),"Closing netcdf file.");
  if (debug>0) fprintf(stdout,"Successfully closed ncid %i\n",ncid);
  ncid=-1;
  frameID=-1; ncframe=-1; atomID=-1; ncatom=-1; coordID=-1;
  cellLengthID=-1; cellAngleID=-1;
  return;
}

/*
 * AmberNetcdf::open()
 * Open up Netcdf file and set all dimension and variable IDs.
 * This is done every time the file is opened up since Im not sure
 * the variable IDs stay the same throughout each opening.
 * Could eventually be separated. 
 */
int AmberNetcdf::open() {
  //fprintf(stdout,"DEBUG: AmberNetcdf::open() called for %s, ncid=%i\n",File->filename,ncid);
  // If already open, return
  if (ncid!=-1) return 0;

  switch (File->access) {

    case READ :
      if (checkNCerr(nc_open(File->filename,NC_NOWRITE,&ncid),
          "Opening Netcdf file %s for reading",File->filename)!=0) return 1;
      break;
    
    case APPEND: 
    case WRITE:
      if (checkNCerr(nc_open(File->filename,NC_WRITE,&ncid),
          "Opening Netcdf file %s for Write.\n",File->filename)!=0) return 1;

  }

  if (debug>0) fprintf(stdout,"Successfully opened %s, ncid=%i\n",File->filename,ncid);
  if (debug>1) this->Debug();
  // Netcdf files are always seekable
  seekable=1;

  frameID=GetDimInfo(NCFRAME,&ncframe);
  if (frameID==-1) return 1;
  atomID=GetDimInfo(NCATOM,&ncatom);
  if (atomID==-1) return 1;
  if (checkNCerr(nc_inq_varid(ncid,NCCOORDS,&coordID),
      "Getting coordinate ID")!=0) return 1;

  // Box info 
  if ( nc_inq_varid(ncid,"cell_lengths",&cellLengthID)==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(ncid,"cell_angles",&cellAngleID),
      "Getting cell angles.")!=0) return 1;
    if (debug>0) fprintf(stdout,"  Netcdf Box information found.\n"); 
    if (P->ifbox==0) {
      fprintf(stderr,"Warning: Netcdf file contains box info but no box info found\n");
      fprintf(stderr,"         in associated parmfile %s; defaulting to rectangular\n",
              P->parmName);
      fprintf(stderr,"         box.\n");
      isBox=1;
    } else {
      isBox=P->ifbox;
    }
  } 

  // Replica Temperatures
  if ( nc_inq_varid(ncid,"temp0",&TempVarID) == NC_NOERR ) {
    if (debug>0) fprintf(stdout,"    Netcdf file has replica temperatures.\n");
    hasTemperature=1;
  } else 
    TempVarID=-1;

  // NOTE: TO BE ADDED
  //spatialID, labelID;
  //int cell_spatialID, cell_angularID;
  //int spatialVID, timeVID, cell_spatialVID, cell_angularVID;

  return 0;
}

/*
 * AmberNetcdf::SetupRead()
 * Just a frontend to open for now. Also check number of atoms.
 */
int AmberNetcdf::SetupRead() {
  if (open()) return 1;
  if (ncatom!=P->natom) {
    fprintf(stdout,"Warning: Number of atoms in NetCDF file %s (%i) does not\n",
            File->filename,ncatom);
    fprintf(stdout,"         match those in associated parmtop (%i)!\n",P->natom);
    return 1;
  }
  stop=ncframe;
  Frames=stop; 
  Coord=(float*) malloc(ncatom*3*sizeof(float));
  close();
  return 0;
}

/*
 * AmberNetcdf::SetupWrite()
 * Create Netcdf file specified by filename and set it up. 
 */
int AmberNetcdf::SetupWrite() {
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[3], count[3];
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };

  // Create file
  if (checkNCerr(nc_create(File->filename,NC_64BIT_OFFSET,&ncid),
    "Creating Netcdf file %s",File->filename)) return 1;
  if (debug>0) 
    fprintf(stdout,"    Successfully created Netcdf file %s, ncid %i\n",File->filename,ncid);

  // Frame, Time
  if (checkNCerr(nc_def_dim(ncid,NCFRAME,NC_UNLIMITED,&frameID),
    "Defining frame dimension.")) return 1;
  dimensionID[0]=frameID;
  if (checkNCerr(nc_def_var(ncid,NCTIME,NC_FLOAT,1,dimensionID,&timeVID),
    "Defining time variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,timeVID,"units",10,"picosecond"),
    "Writing time VID units.")) return 1;

  // Spatial
  if (checkNCerr(nc_def_dim(ncid,NCSPATIAL,3,&spatialID),
    "Defining spatial dimension.")) return 1;
  dimensionID[0] = spatialID;
  if (checkNCerr(nc_def_var(ncid,NCSPATIAL,NC_CHAR,1,dimensionID,&spatialVID),
    "Defining spatial variable.")) return 1;

  // Atoms
  if (checkNCerr(nc_def_dim(ncid,NCATOM,P->natom,&atomID),
    "Defining atom dimension.")) return 1;
  ncatom=P->natom;

  // Coords
  dimensionID[0] = frameID;
  dimensionID[1] = atomID;
  dimensionID[2] = spatialID;
  if (checkNCerr(nc_def_var(ncid,NCCOORDS,NC_FLOAT,3,dimensionID,&coordID),
    "Defining coordinates variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,coordID,"units",8,"angstrom"),
    "Writing coordinates variable units.")) return 1;

  // Cell Spatial
  if (checkNCerr(nc_def_dim(ncid,NCCELL_SPATIAL,3,&cell_spatialID),
    "Defining cell spatial dimension.")) return 1;
  dimensionID[0]=cell_spatialID;
  if (checkNCerr(nc_def_var(ncid,NCCELL_SPATIAL,NC_CHAR,1,dimensionID,&cell_spatialVID),
    "Defining cell spatial variable.")) return 1;

  // Cell angular
  if (checkNCerr(nc_def_dim(ncid,NCLABEL,NCLABELLEN,&labelID),
    "Defining label dimension.")) return 1;
  if (checkNCerr(nc_def_dim(ncid,NCCELL_ANGULAR,3,&cell_angularID),
    "Defining cell angular dimension.")) return 1;
  dimensionID[0] = cell_angularID;
  dimensionID[1] = labelID;
  if (checkNCerr(nc_def_var(ncid,NCCELL_ANGULAR,NC_CHAR,2,dimensionID,&cell_angularVID),
    "Defining cell angular variable.")) return 1;

  // Box Info
  if (isBox>0) {
    dimensionID[0]=frameID;
    dimensionID[1]=cell_spatialID;
    if (checkNCerr(nc_def_var(ncid,"cell_lengths",NC_DOUBLE,2,dimensionID,&cellLengthID),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellLengthID,"units",8,"angstrom"),
    "Writing cell length variable units.")) return 1;
    dimensionID[1]=cell_angularID;
    if (checkNCerr(nc_def_var(ncid,"cell_angles",NC_DOUBLE,2,dimensionID,&cellAngleID),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellAngleID,"units",6,"degree"),
    "Writing cell angle variable units.")) return 1;
  }

  // Set up title
  if (title==NULL) {
    title=(char*) malloc(30*sizeof(char));
    strcpy(title,"Cpptraj Generated trajectory");
  }

  // Attributes
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"title",strlen(title),title),
    "Writing title.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"application",5,"AMBER"),
    "Writing application.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"program",7,"cpptraj"),
    "Writing program.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"programVersion",3,"1.0"),
    "Writing program version.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"Conventions",5,"AMBER"),
    "Writing conventions.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"ConventionVersion",3,"1.0"),
    "Writing conventions version.")) return 1;

  /* DAN ROE: Replica temperature 
  if (trajInfo->isREMDTRAJ) {
    fprintf(stdout,"NETCDF: Defining replica temperature in output trajectory.\n");
    dimensionID[0] = NCInfo->frameDID;
    netcdfDefineVariable(NCInfo->ncid, "temp0", NC_DOUBLE, 1, dimensionID, &NCInfo->TempVarID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->TempVarID,"units","kelvin");
  }*/

  // Set fill mode
  if (checkNCerr(nc_set_fill(ncid, NC_NOFILL, dimensionID),
    "NetCDF setting fill value.")) return 1;

  // End netcdf definitions
  if (checkNCerr(nc_enddef(ncid),"NetCDF error on ending definitions."))
    return 1;

  // Specify spatial dimension labels
  start[0] = 0;
  count[0] = 3;

  xyz[0] = 'x'; xyz[1] = 'y'; xyz[2] = 'z';
  if (checkNCerr(nc_put_vara_text(ncid, spatialVID, start, count, xyz),
    "Error on NetCDF output of spatial VID 'x', 'y' and 'z'")) return 1;

  xyz[0] = 'a'; xyz[1] = 'b'; xyz[2] = 'c';
  if (checkNCerr(nc_put_vara_text(ncid, cell_spatialVID, start, count, xyz),
    "Error on NetCDF output of cell spatial VID 'a', 'b' and 'c'")) return 1;

  start[0] = 0; start[1] = 0;
  count[0] = 3; count[1] = 5;
  if (checkNCerr(nc_put_vara_text(ncid, cell_angularVID, start, count, abc),
    "Error on NetCDF output of cell angular VID 'alpha', 'beta ' and 'gamma'")) 
    return 1;

  // Allocate memory and close the file. It will be reopened WRITE
  Coord=(float*) malloc(ncatom*3*sizeof(float));
  close();
  
  return 0;
}


/*
 * AmberNetcdf::getFrame()
 * Get the specified frame from amber netcdf file
 * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
 */
int AmberNetcdf::getFrame(int set) {
  size_t start[3], count[3];

  // Get temperature
  if (TempVarID!=-1) {
    start[0]=set;
    count[0]=1;
    if ( checkNCerr(nc_get_vara_double(ncid, TempVarID, start,count,&(F->T)),
                    "Getting replica temperature.")!=0 ) return 1;
    //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
  }

  // Read Coords 
  start[0]=set;
  start[1]=0;
  start[2]=0;
  count[0]=1;
  count[1]=ncatom;
  count[2]=3;
  if ( checkNCerr(nc_get_vara_float(ncid, coordID, start, count, Coord),
                  "Getting frame %i",set)!=0 ) return 1;

  // Read box info 
  if (isBox!=0) {
    count [1]=3;
    count [2]=0;
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthID, start, count, F->box),
                    "Getting cell lengths.")!=0 ) return 1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleID, start, count, &(F->box[3])),
                    "Getting cell angles.")!=0 ) return 1;
  }

  F->floatToFrame(Coord);

  return 0;
}

/*
 * AmberNetcdf::writeFrame() 
 */
int AmberNetcdf::writeFrame(int set) {
  size_t start[3], count[3];

  // If coords have been stripped and P->natom is different from ncatom, need
  // to redefine atom IDs etc.
  // NOTE: No checking is done to see if this happens multiple times, which
  //       will really screw things up.

  //if (NC==NULL) {
  //  fprintf(stderr,"Error in writeNCframe: output file not open.");
  //  return 1;
  //}
  F->frameToFloat(Coord);

  // write coords
  start[0]=set;
  start[1]=0;
  start[2]=0;
  count[0]=1;
  count[1]=F->natom; // Use F->natom in case of strip
  count[2]=3;
  if (checkNCerr(nc_put_vara_float(ncid,coordID,start,count,Coord),
    "Netcdf Writing frame %i",set)) return 1;

  // write box
  if (isBox>0 && cellLengthID!=0) {
    count[1]=3;
    count[2]=0;
    if (checkNCerr(nc_put_vara_double(ncid,cellLengthID,start,count,F->box),
      "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(ncid,cellAngleID,start,count, &(F->box[3])),
      "Writing cell angles.")) return 1;
  }
  
  nc_sync(ncid); // Necessary after every write??

  return 0;
}  

/*
 * Info()
 */
void AmberNetcdf::Info() {
  fprintf(stdout,"  File (%s) is a NetCDF AMBER trajectory", File->filename
            //(p->isVelocity ? " and velocities" : "")
         );

  if (TempVarID!=-1) fprintf(stdout," with replica temperatures");

  /*if (debug > 2) {
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
