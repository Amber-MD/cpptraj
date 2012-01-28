#ifdef BINTRAJ
// This file contains a collection of routines designed for reading
// netcdf trajectory files used with amber.
// Dan Roe 10-2008
// Original implementation of netcdf in Amber by Jon Mongan.
#include <cstdlib>
#include <cstring> // For title length
#include "netcdf.h"
#include "Traj_AmberNetcdf.h"
#include "NetcdfRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AmberNetcdf::AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Constructor\n");
  ncid=-1;
  frameDID=-1;
  ncframe=-1;
  atomDID=-1;
  ncatom=-1;
  ncatom3=-1;
  Coord=NULL;
  coordVID=-1;
  cellAngleVID=-1;
  cellLengthVID=-1;

  spatialDID=-1;
  labelDID=-1;
  cell_spatialDID=-1;
  cell_angularDID=-1;
  spatialVID=-1;
  timeVID=-1;
  cell_spatialVID=-1;
  cell_angularVID=-1;
  TempVID=-1;

  remd_dimension=0;
  dimensionDID=-1;
  groupnumVID=-1;
  dimtypeVID=-1;
  indicesVID=-1;
  remd_groupnum=NULL;
  remd_dimtype=NULL;
  remd_indices=NULL;

  // Netcdf files are always seekable
  seekable=true;
} 

// DESTRUCTOR
AmberNetcdf::~AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  this->closeTraj();
  if (Coord!=NULL) free(Coord);
  if (remd_groupnum!=NULL) delete[] remd_groupnum;
  if (remd_dimtype!=NULL) delete[] remd_dimtype;
  if (remd_indices!=NULL) delete[] remd_indices;
  // NOTE: Need to close file?
}

// AmberNetcdf::close()
/** Close netcdf file. Set ncid to -1 since it can change between open
  * and close calls.
  */
void AmberNetcdf::closeTraj() {
  if (ncid<0) return;
  //if (tfile->access!=READ) nc_sync(ncid); 
  checkNCerr(nc_close(ncid),"Closing netcdf file.");
  if (debug>0) rprintf("Successfully closed ncid %i\n",ncid);
  ncid=-1;
  return;
}

// AmberNetcdf::openTraj()
/** Open up Netcdf file and set ncid. Variable and Dimension IDs are set up
  * by SetupRead / SetupWrite and will not change for a given file between
  * open and close calls.
  */
int AmberNetcdf::openTraj() {
  //fprintf(stdout,"DEBUG: AmberNetcdf::openTraj() called for %s, ncid=%i\n",tfile->filename,ncid);
  // If already open, return
  if (ncid!=-1) return 0;

  switch (tfile->access) {
    case READ :
      if (checkNCerr(nc_open(tfile->filename,NC_NOWRITE,&ncid),
          "Opening Netcdf file %s for reading",tfile->filename)!=0) return 1;
      break;
    
    case APPEND: 
    case WRITE:
      if (checkNCerr(nc_open(tfile->filename,NC_WRITE,&ncid),
          "Opening Netcdf file %s for Write.\n",tfile->filename)!=0) return 1;
  }

  if (debug>0) rprintf("Successfully opened %s, ncid=%i\n",tfile->filename,ncid);
  if (debug>1) NetcdfDebug(ncid);

  return 0;
}

// AmberNetcdf::setupRead()
/** Open the netcdf file, read all dimension and variable IDs, close.
  * Return the number of frames in the file. 
  */
int AmberNetcdf::setupRead(AmberParm* trajParm) {
  char *attrText;            // For checking conventions and version 
  int spatial;               // For checking spatial dimensions
  double box[6];             // For checking box type
  size_t start[3], count[3]; // For checking box type

  if (openTraj()) return -1;

  // Get global attributes
  if (title==NULL) title = GetAttrText(ncid,NC_GLOBAL, "title");
  attrText = GetAttrText(ncid,NC_GLOBAL, "Conventions");
  if (attrText==NULL || strstr(attrText,"AMBER")==NULL) 
    rprintf("WARNING: Netcdf file %s conventions do not include \"AMBER\" (%s)\n",
            tfile->filename, attrText);
  if (attrText!=NULL) free(attrText);
  attrText = GetAttrText(ncid,NC_GLOBAL, "ConventionVersion");
  if (attrText==NULL || strcmp(attrText,"1.0")!=0)
    rprintf("WARNING: Netcdf file %s has ConventionVersion that is not 1.0 (%s)\n",
            tfile->filename, attrText);
  if (attrText!=NULL) free(attrText);

  // Get frame, atoms, coord, and spatial info
  frameDID=GetDimInfo(ncid,NCFRAME,&ncframe);
  if (frameDID==-1) return -1;
  atomDID=GetDimInfo(ncid,NCATOM,&ncatom);
  if (atomDID==-1) return -1;
  ncatom3 = ncatom * 3;
  if (checkNCerr(nc_inq_varid(ncid,NCCOORDS,&coordVID),
      "Getting coordinate ID")!=0) return -1;
  attrText = GetAttrText(ncid,coordVID, "units");
  if (attrText==NULL || strcmp(attrText,"angstrom")!=0) 
    rprintf("WARNING: Netcdf file %s has length units of %s - expected angstrom.\n",
            tfile->filename,attrText);
  if (attrText!=NULL) free(attrText);
  spatialDID=GetDimInfo(ncid,NCSPATIAL,&spatial);
  if (spatialDID==-1) return -1;
  if (spatial!=3) {
    rprintf("Error: ncOpen: Expected 3 spatial dimenions in %s, got %i\n",
            tfile->filename, spatial);
    return -1;
  }
  if ( checkNCerr(nc_inq_varid(ncid, NCSPATIAL, &spatialVID),
       "Getting spatial VID\n")) return -1;

  // Sanity check on Time
  if ( checkNCerr( nc_inq_varid(ncid, NCTIME, &timeVID),
       "Getting Netcdf time VID.")) return -1;
  attrText = GetAttrText(ncid,timeVID, "units");
  if (attrText==NULL || strcmp(attrText,"picosecond")!=0) 
    rprintf("WARNING: Netcdf file %s has time units of %s - expected picosecond.\n",
            tfile->filename, attrText);
  if (attrText!=NULL) free(attrText);

  // Box info
  // NOTE: If no box info found in parm should really try to determine correct
  //       box type from angles. 
  if ( nc_inq_varid(ncid,NCCELL_LENGTHS,&cellLengthVID)==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(ncid,NCCELL_ANGLES,&cellAngleVID),
      "Getting cell angles.")!=0) return -1;
    if (debug>0) mprintf("  Netcdf Box information found.\n");
    hasBox=true;
    // If present, get box angles. These will be used to determine the box 
    // type in TrajectoryFile.
    start[0]=0; start[1]=0; start[2]=0; 
    count[0]=1; count[1]=3; count[2]=0;
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthVID, start, count, box),
                    "Getting cell lengths.")!=0 ) return -1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleVID, start, count, boxAngle),
                    "Getting cell angles.")!=0 ) return -1;
  } 

  // Replica Temperatures
  if ( nc_inq_varid(ncid,NCTEMPERATURE,&TempVID) == NC_NOERR ) {
    if (debug>0) mprintf("    Netcdf file has replica temperatures.\n");
    hasTemperature=true;
  } else 
    TempVID=-1;

  // Multi-d REMD info
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
  }

  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;

  // Check that specified number of atoms matches expected number.
  if (ncatom!=trajParm->natom) {
    mprinterr("Error: Number of atoms in NetCDF file %s (%i) does not\n",
              tfile->filename,ncatom);
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->natom);
    return -1;
  }

  // Amber Netcdf coords are float. Allocate a float array for converting
  // float to/from double.
  Coord=(float*) realloc(Coord, ncatom3*sizeof(float));
  closeTraj();
  return ncframe;
}

// AmberNetcdf::processWriteArgs()
int AmberNetcdf::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("remdtraj")) this->SetTemperature();
  return 0;
}

// AmberNetcdf::setupWrite()
/** Create Netcdf file specified by filename and set up dimension and
  * variable IDs. 
  */
int AmberNetcdf::setupWrite(AmberParm *trajParm) {
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[3], count[3];
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };

  // Create file
  if (checkNCerr(nc_create(tfile->filename,NC_64BIT_OFFSET,&ncid),
    "Creating Netcdf file %s",tfile->filename)) return 1;
  if (debug>0) 
    mprintf("    Successfully created Netcdf file %s, ncid %i\n",tfile->filename,ncid);

  // Frame, Time
  ncframe=0;
  if (checkNCerr(nc_def_dim(ncid,NCFRAME,NC_UNLIMITED,&frameDID),
    "Defining frame dimension.")) return 1;
  dimensionID[0]=frameDID;
  if (checkNCerr(nc_def_var(ncid,NCTIME,NC_FLOAT,1,dimensionID,&timeVID),
    "Defining time variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,timeVID,"units",10,"picosecond"),
    "Writing time VID units.")) return 1;

  // Spatial
  if (checkNCerr(nc_def_dim(ncid,NCSPATIAL,3,&spatialDID),
    "Defining spatial dimension.")) return 1;
  dimensionID[0] = spatialDID;
  if (checkNCerr(nc_def_var(ncid,NCSPATIAL,NC_CHAR,1,dimensionID,&spatialVID),
    "Defining spatial variable.")) return 1;

  // Atoms
  if (checkNCerr(nc_def_dim(ncid,NCATOM,trajParm->natom,&atomDID),
    "Defining atom dimension.")) return 1;
  ncatom=trajParm->natom;
  ncatom3 = ncatom * 3;

  // Coords
  dimensionID[0] = frameDID;
  dimensionID[1] = atomDID;
  dimensionID[2] = spatialDID;
  if (checkNCerr(nc_def_var(ncid,NCCOORDS,NC_FLOAT,3,dimensionID,&coordVID),
    "Defining coordinates variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,coordVID,"units",8,"angstrom"),
    "Writing coordinates variable units.")) return 1;

  // Cell Spatial
  if (checkNCerr(nc_def_dim(ncid,NCCELL_SPATIAL,3,&cell_spatialDID),
    "Defining cell spatial dimension.")) return 1;
  dimensionID[0]=cell_spatialDID;
  if (checkNCerr(nc_def_var(ncid,NCCELL_SPATIAL,NC_CHAR,1,dimensionID,&cell_spatialVID),
    "Defining cell spatial variable.")) return 1;

  // Cell angular
  // NOTE: Should this even be defined if nobox?
  if (checkNCerr(nc_def_dim(ncid,NCLABEL,NCLABELLEN,&labelDID),
    "Defining label dimension.")) return 1;
  if (checkNCerr(nc_def_dim(ncid,NCCELL_ANGULAR,3,&cell_angularDID),
    "Defining cell angular dimension.")) return 1;
  dimensionID[0] = cell_angularDID;
  dimensionID[1] = labelDID;
  if (checkNCerr(nc_def_var(ncid,NCCELL_ANGULAR,NC_CHAR,2,dimensionID,&cell_angularVID),
    "Defining cell angular variable.")) return 1;

  // Box Info
  if (hasBox) {
    dimensionID[0]=frameDID;
    dimensionID[1]=cell_spatialDID;
    if (checkNCerr(nc_def_var(ncid,NCCELL_LENGTHS,NC_DOUBLE,2,dimensionID,&cellLengthVID),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellLengthVID,"units",8,"angstrom"),
    "Writing cell length variable units.")) return 1;
    dimensionID[1]=cell_angularDID;
    if (checkNCerr(nc_def_var(ncid,NCCELL_ANGLES,NC_DOUBLE,2,dimensionID,&cellAngleVID),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellAngleVID,"units",6,"degree"),
    "Writing cell angle variable units.")) return 1;
  }

  // Set up title
  if (title==NULL)
    this->SetTitle((char*)"Cpptraj Generated trajectory\0");

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

  // Replica temperature
  if (hasTemperature) { 
    mprintf("NETCDF: Defining replica temperature in output trajectory.\n");
    dimensionID[0] = frameDID;
    if (checkNCerr(nc_def_var(ncid,NCTEMPERATURE,NC_DOUBLE,1,dimensionID,&TempVID),
        "NetCDF error on defining replica temperature")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,TempVID,"units",6,"kelvin"),
        "NetCDF error on defining temperature units.")) return 1;
  }

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
  Coord=(float*) realloc(Coord, ncatom3*sizeof(float));
  closeTraj();
  
  return 0;
}

// AmberNetcdf::readFrame()
/** Get the specified frame from amber netcdf file
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  */
int AmberNetcdf::readFrame(int set,double *X, double *V,double *box, double *T) {
  size_t start[3], count[3];

  // Get temperature
  if (TempVID!=-1) {
    start[0]=set;
    count[0]=1;
    if ( checkNCerr(nc_get_vara_double(ncid, TempVID, start,count,T),
                    "Getting replica temperature.")!=0 ) return 1;
    //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",F->T);
  }

  // Get replica indices
  if (indicesVID!=-1) {
    start[0]=set;
    start[1]=0;
    count[0]=1;
    count[1]=remd_dimension;
    if ( checkNCerr(nc_get_vara_int(ncid, indicesVID, start, count, remd_indices),
                    " Getting replica indices.")!=0 ) return 1;
    mprintf("DEBUG:\tReplica indices:");
    for (int dim=0; dim < remd_dimension; dim++) mprintf(" %i",remd_indices[dim]);
    mprintf("\n");
  }

  // Read Coords 
  start[0]=set;
  start[1]=0;
  start[2]=0;
  count[0]=1;
  count[1]=ncatom;
  count[2]=3;
  if ( checkNCerr(nc_get_vara_float(ncid, coordVID, start, count, Coord),
                  "Getting frame %i",set)!=0 ) return 1;

  // Read box info 
  if (hasBox) {
    count [1]=3;
    count [2]=0;
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthVID, start, count, box),
                    "Getting cell lengths.")!=0 ) return 1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleVID, start, count, box+3),
                    "Getting cell angles.")!=0 ) return 1;
  }

  FloatToDouble(X,Coord,ncatom3);

  return 0;
}

// AmberNetcdf::writeFrame() 
int AmberNetcdf::writeFrame(int set, double *X, double *V, double *box, double T) {
  size_t start[3], count[3];

  DoubleToFloat(Coord,X,ncatom3);

  // Write coords
  start[0]=ncframe;
  start[1]=0;
  start[2]=0;
  count[0]=1;
  count[1]=ncatom;
  count[2]=3;
  if (checkNCerr(nc_put_vara_float(ncid,coordVID,start,count,Coord),
    "Netcdf Writing frame %i",set)) return 1;

  // Write box
  if (hasBox) {
    count[1]=3;
    count[2]=0;
    if (checkNCerr(nc_put_vara_double(ncid,cellLengthVID,start,count,box),
      "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(ncid,cellAngleVID,start,count, box+3),
      "Writing cell angles.")) return 1;
  }

  // Write temperature
  if (TempVID!=-1) {
    if ( checkNCerr( nc_put_vara_double(ncid,TempVID,start,count,&T),
         "Writing temperature.") ) return 1;
  }
  
  nc_sync(ncid); // Necessary after every write??

  ncframe++;

  return 0;
}  

// AmberNetcdf::info()
void AmberNetcdf::info() {
  mprintf("is a NetCDF AMBER trajectory"
            //(p->isVelocity ? " and velocities" : "")
         );

  if (hasTemperature) mprintf(" with replica temperatures");

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
