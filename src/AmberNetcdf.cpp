#ifdef BINTRAJ
/* 
 * This file contains a collection of routines designed for reading
 * (and writing?) netcdf trajectory files used with amber.
 * Dan Roe 10-2008
 * Original implementation of netcdf in Amber by Jon Mongan.
 */
#include <cstdlib>
#include <cstring> // For title length
#include "netcdf.h"
#include "AmberNetcdf.h"
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
} 

// DESTRUCTOR
AmberNetcdf::~AmberNetcdf() {
  //fprintf(stderr,"Amber Netcdf Destructor\n");
  if (Coord!=NULL) free(Coord);
  // NOTE: Need to close file?
}

/*
 * AmberNetcdf::close()
 */
void AmberNetcdf::close() {
  checkNCerr(nc_close(ncid),"Closing netcdf file.");
  if (debug>0) rprintf("Successfully closed ncid %i\n",ncid);
  ncid=-1;
  frameDID=-1; ncframe=-1; atomDID=-1; ncatom=-1; coordVID=-1;
  cellLengthVID=-1; cellAngleVID=-1;
  return;
}

/*
 * AmberNetcdf::open()
 * Open up Netcdf file and set all dimension and variable IDs.
 * This is done every time the file is opened up since Im not sure
 * the variable IDs stay the same throughout each opening.
 * Could eventually be separated.
 * NOTE: Replace attrText allocs with static buffer? 
 */
int AmberNetcdf::open() {
  char *attrText; // For checking conventions and version 
  int spatial; // For checking spatial dimensions

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

  if (debug>0) rprintf("Successfully opened %s, ncid=%i\n",File->filename,ncid);
  if (debug>1) NetcdfDebug(ncid);
  // Netcdf files are always seekable
  seekable=1;

  // Get global attributes
  if (title==NULL) title = GetAttrText(ncid,NC_GLOBAL, "title");
  attrText = GetAttrText(ncid,NC_GLOBAL, "Conventions");
  if (attrText==NULL || strstr(attrText,"AMBER")==NULL) 
    rprintf("WARNING: Netcdf file %s conventions do not include \"AMBER\" (%s)\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);
  attrText = GetAttrText(ncid,NC_GLOBAL, "ConventionVersion");
  if (attrText==NULL || strcmp(attrText,"1.0")!=0)
    rprintf("WARNING: Netcdf file %s has ConventionVersion that is not 1.0 (%s)\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);

  // Get frame, atoms, coord, and spatial info
  frameDID=GetDimInfo(ncid,NCFRAME,&ncframe);
  if (frameDID==-1) return 1;
  atomDID=GetDimInfo(ncid,NCATOM,&ncatom);
  if (atomDID==-1) return 1;
  if (checkNCerr(nc_inq_varid(ncid,NCCOORDS,&coordVID),
      "Getting coordinate ID")!=0) return 1;
  attrText = GetAttrText(ncid,coordVID, "units");
  if (attrText==NULL || strcmp(attrText,"angstrom")!=0) 
    rprintf("WARNING: Netcdf file %s has length units of %s - expected angstrom.\n",
            File->filename,attrText);
  if (attrText!=NULL) free(attrText);
  spatialDID=GetDimInfo(ncid,NCSPATIAL,&spatial);
  if (spatialDID==-1) return 1;
  if (spatial!=3) {
    rprintf("Error: ncOpen: Expected 3 spatial dimenions in %s, got %i\n",
            File->filename, spatial);
    return 1;
  }
  if ( checkNCerr(nc_inq_varid(ncid, NCSPATIAL, &spatialVID),
       "Getting spatial VID\n")) return 1;

  // Sanity check on Time
  if ( checkNCerr( nc_inq_varid(ncid, NCTIME, &timeVID),
       "Getting Netcdf time VID.")) return 1;
  attrText = GetAttrText(ncid,timeVID, "units");
  if (attrText==NULL || strcmp(attrText,"picosecond")!=0) 
    rprintf("WARNING: Netcdf file %s has time units of %s - expected picosecond.\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);

  // Box info
  // NOTE: If no box info found in parm should really try to determine correct
  //       box type from angles. 
  if ( nc_inq_varid(ncid,"cell_lengths",&cellLengthVID)==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(ncid,"cell_angles",&cellAngleVID),
      "Getting cell angles.")!=0) return 1;
    if (debug>0) mprintf("  Netcdf Box information found.\n"); 
    if (P->ifbox==0) {
      mprintf("Warning: Netcdf file contains box info but no box info found\n");
      mprintf("         in associated parmfile %s; defaulting to orthogonal.\n",
              P->parmName);
      isBox=1;
    } else {
      isBox=P->ifbox;
    }
  } 

  // Replica Temperatures
  if ( nc_inq_varid(ncid,NCTEMPERATURE,&TempVID) == NC_NOERR ) {
    if (debug>0) mprintf("    Netcdf file has replica temperatures.\n");
    hasTemperature=1;
  } else 
    TempVID=-1;

  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;

  return 0;
}

/*
 * AmberNetcdf::SetupRead()
 * Just a frontend to open for now. Also check number of atoms.
 */
int AmberNetcdf::SetupRead() {
  if (open()) return 1;
  if (ncatom!=P->natom) {
    rprintf("Warning: Number of atoms in NetCDF file %s (%i) does not\n",
            File->filename,ncatom);
    rprintf("         match those in associated parmtop (%i)!\n",P->natom);
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
    mprintf("    Successfully created Netcdf file %s, ncid %i\n",File->filename,ncid);

  // Frame, Time
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
  if (checkNCerr(nc_def_dim(ncid,NCATOM,P->natom,&atomDID),
    "Defining atom dimension.")) return 1;
  ncatom=P->natom;

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
  if (checkNCerr(nc_def_dim(ncid,NCLABEL,NCLABELLEN,&labelDID),
    "Defining label dimension.")) return 1;
  if (checkNCerr(nc_def_dim(ncid,NCCELL_ANGULAR,3,&cell_angularDID),
    "Defining cell angular dimension.")) return 1;
  dimensionID[0] = cell_angularDID;
  dimensionID[1] = labelDID;
  if (checkNCerr(nc_def_var(ncid,NCCELL_ANGULAR,NC_CHAR,2,dimensionID,&cell_angularVID),
    "Defining cell angular variable.")) return 1;

  // Box Info
  if (isBox>0) {
    dimensionID[0]=frameDID;
    dimensionID[1]=cell_spatialDID;
    if (checkNCerr(nc_def_var(ncid,"cell_lengths",NC_DOUBLE,2,dimensionID,&cellLengthVID),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellLengthVID,"units",8,"angstrom"),
    "Writing cell length variable units.")) return 1;
    dimensionID[1]=cell_angularDID;
    if (checkNCerr(nc_def_var(ncid,"cell_angles",NC_DOUBLE,2,dimensionID,&cellAngleVID),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellAngleVID,"units",6,"degree"),
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
    mprintf("NETCDF: Defining replica temperature in output trajectory.\n");
    dimensionID[0] = NCInfo->frameDID;
    netcdfDefineVariable(NCInfo->ncid, NCTEMPERATURE, NC_DOUBLE, 1, dimensionID, &NCInfo->TempVID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->TempVID,"units","kelvin");
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
  if (TempVID!=-1) {
    start[0]=set;
    count[0]=1;
    if ( checkNCerr(nc_get_vara_double(ncid, TempVID, start,count,&(F->T)),
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
  if ( checkNCerr(nc_get_vara_float(ncid, coordVID, start, count, Coord),
                  "Getting frame %i",set)!=0 ) return 1;

  // Read box info 
  if (isBox!=0) {
    count [1]=3;
    count [2]=0;
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthVID, start, count, F->box),
                    "Getting cell lengths.")!=0 ) return 1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleVID, start, count, &(F->box[3])),
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
  // NOTE: This should not occur since all actions are completed before the write
  //       stage. However if a write action is ever implemented this could become
  //       a problem.
  //if (P->natom != ncatom) {
  //  fprintf(stdout,"Error: AmberNetcdf::writeFrame(%i)\n",set);
  //  fprintf(stdout,"       %s was set up for %i atoms but attempting to write %i atoms!\n",
  //          trajfilename, ncatom, P->natom);
  //  fprintf(stdout,"       This can happen e.g. if there are multiple strip commands\n");
  //  fprintf(stdout,"       and is not supported by the Netcdf file format.\n");
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
  if (checkNCerr(nc_put_vara_float(ncid,coordVID,start,count,Coord),
    "Netcdf Writing frame %i",set)) return 1;

  // write box
  if (isBox>0 && cellLengthVID!=-1) {
    count[1]=3;
    count[2]=0;
    if (checkNCerr(nc_put_vara_double(ncid,cellLengthVID,start,count,F->box),
      "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(ncid,cellAngleVID,start,count, &(F->box[3])),
      "Writing cell angles.")) return 1;
  }
  
  nc_sync(ncid); // Necessary after every write??

  return 0;
}  

/*
 * Info()
 */
void AmberNetcdf::Info() {
  mprintf("  File (%s) is a NetCDF AMBER trajectory", File->filename
            //(p->isVelocity ? " and velocities" : "")
         );

  if (TempVID!=-1) mprintf(" with replica temperatures");

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
