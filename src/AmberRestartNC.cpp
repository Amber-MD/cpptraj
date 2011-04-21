#ifdef BINTRAJ
/*
 * EXPERIMENTAL 
 * This file contains a collection of routines designed for reading
 * (and writing?) netcdf restart files used with amber.
 * Dan Roe 2011-01-07
 */
#include <cstdio> // sprintf
#include <cstdlib>
#include <cstring> // For title length
#include "netcdf.h"
#include "AmberRestartNC.h"
#include "NetcdfRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AmberRestartNC::AmberRestartNC() {
  //fprintf(stderr,"Amber Netcdf Restart Constructor\n");
  ncid=-1;
  atomDID=-1;
  ncatom=-1;
  coordVID=-1;
  hasVelocity=0;
  velocityVID=-1;
  velocityScale=20.455; // Will be overwritten when read in but should be the same 
  cellAngleVID=-1;
  cellLengthVID=-1;

  spatialDID=-1;
  labelDID=-1;
  cell_spatialDID=-1;
  cell_angularDID=-1;
  spatialVID=-1;
  timeVID=-1;
  restartTime=0.0;
  cell_spatialVID=-1;
  cell_angularVID=-1;
  TempVID=-1;
} 

// DESTRUCTOR
AmberRestartNC::~AmberRestartNC() {
  //fprintf(stderr,"Amber Netcdf Restart Destructor\n");
  // NOTE: Need to close file?
}

/*
 * AmberRestartNC::close()
 */
void AmberRestartNC::close() {
  if (ncid==-1) return;
  checkNCerr(nc_close(ncid),"Closing netcdf Restart file.");
  if (debug>0) mprintf("Successfully closed ncid %i\n",ncid);
  ncid=-1;
  atomDID=-1; ncatom=-1; coordVID=-1; velocityVID=-1;
  cellLengthVID=-1; cellAngleVID=-1;
  return;
}

/*
 * AmberRestartNC::open()
 * Open up Netcdf restart file and set all dimension and variable IDs.
 */
int AmberRestartNC::open() {
  mprintf("DEBUG: AmberRestartNC::open() called for %s, ncid=%i\n",File->filename,ncid);
  // If already open, return
  if (ncid!=-1) return 0;

  switch (File->access) {

    case READ :
      if (checkNCerr(nc_open(File->filename,NC_NOWRITE,&ncid),
          "Opening Netcdf restart file %s for reading",File->filename)!=0) return 1;
      break;
    
    case APPEND: 
      mprintf("Error: %s - Append is not supported by netcdf restart.\n",File->filename);
      return 1;
      break;
    case WRITE:
      // All opening and closing of writes are done per-frame
      return 0;
  }

  if (debug>0) mprintf("Successfully opened restart %s, ncid=%i\n",File->filename,ncid);
  if (debug>1) NetcdfDebug(ncid);
  // Netcdf files are always seekable
  seekable=1;

  return 0;
}

/*
 * AmberRestartNC::SetupRead()
 * Set up netcdf restart file for reading, get all variable and dimension IDs. 
 * Also check number of atoms against associated parmtop.
 * NOTE: Replace attrText allocs with static buffer? 
 */
int AmberRestartNC::SetupRead() {
  char *attrText; // For checking conventions and version 
  int spatial; // For checking spatial dimensions
  double box[6];
  size_t start[2], count[2];

  if (open()) return 1;

  // Get global attributes
  if (title==NULL) title = GetAttrText(ncid,NC_GLOBAL, "title");
  attrText = GetAttrText(ncid,NC_GLOBAL, "Conventions");
  if (attrText==NULL || strstr(attrText,"AMBERRESTART")==NULL) 
    mprintf("WARNING: Netcdf restart file %s conventions do not include \"AMBERRESTART\" (%s)\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);
  attrText = GetAttrText(ncid,NC_GLOBAL, "ConventionVersion");
  if (attrText==NULL || strcmp(attrText,"1.0")!=0)
    mprintf("WARNING: Netcdf restart file %s has ConventionVersion that is not 1.0 (%s)\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);

  // Get atoms, coord, and spatial info
  atomDID=GetDimInfo(ncid,NCATOM,&ncatom);
  if (atomDID==-1) return 1;
  if (checkNCerr(nc_inq_varid(ncid,NCCOORDS,&coordVID),
      "Getting coordinate ID")!=0) return 1;
  attrText = GetAttrText(ncid,coordVID, "units");
  if (attrText==NULL || strcmp(attrText,"angstrom")!=0) 
    mprintf("WARNING: Netcdf file %s has length units of %s - expected angstrom.\n",
            File->filename,attrText);
  if (attrText!=NULL) free(attrText);
  spatialDID=GetDimInfo(ncid,NCSPATIAL,&spatial);
  if (spatialDID==-1) return 1;
  if (spatial!=3) {
    mprintf("Error: ncOpen: Expected 3 spatial dimenions in %s, got %i\n",
            File->filename, spatial);
    return 1;
  }
  if ( checkNCerr(nc_inq_varid(ncid, NCSPATIAL, &spatialVID),
       "Getting spatial VID\n")) return 1;

  // Get Velocity info
  if ( nc_inq_varid(ncid,NCVELO,&velocityVID)==NC_NOERR ) {
    if (debug>0) mprintf("    Netcdf restart file has velocities.\n");
    hasVelocity=1;
  } else
    velocityVID=-1;

  // Get Time info
  if ( checkNCerr( nc_inq_varid(ncid, NCTIME, &timeVID),
       "Getting Netcdf time VID.")) return 1;
  attrText = GetAttrText(ncid,timeVID, "units");
  if (attrText==NULL || strcmp(attrText,"picosecond")!=0) 
    mprintf("WARNING: Netcdf restart file %s has time units of %s - expected picosecond.\n",
            File->filename, attrText);
  if (attrText!=NULL) free(attrText);
  if ( checkNCerr(nc_get_var_double(ncid, timeVID, &restartTime),
       "Getting netcdf restart time.")) return 1;
  mprintf("    Netcdf restart time= %lf\n",restartTime);

  // Box info
  if ( nc_inq_varid(ncid,"cell_lengths",&cellLengthVID)==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(ncid,"cell_angles",&cellAngleVID),
      "Getting cell angles.")!=0) return 1;
    if (debug>0) mprintf("  Netcdf restart Box information found.\n");
    // Determine box type from angles
    start[0]=0; start[1]=0;
    count[0]=3; count[1]=0; 
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthVID, start, count, box),
                    "Getting cell lengths.")!=0 ) return 1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleVID, start, count, box+3),
                    "Getting cell angles.")!=0 ) return 1;
    CheckBoxType(box);
  } 

  // Replica Temperatures
  if ( nc_inq_varid(ncid,NCTEMPERATURE,&TempVID) == NC_NOERR ) {
    if (debug>0) mprintf("    Netcdf restart file has replica temperature.\n");
    hasTemperature=1;
  } else 
    TempVID=-1;

  // NOTE: TO BE ADDED
  // labelDID;
  //int cell_spatialDID, cell_angularDID;
  //int spatialVID, cell_spatialVID, cell_angularVID;

  if (ncatom!=P->natom) {
    mprintf("Warning: Number of atoms in NetCDF restart file %s (%i) does not\n",
            File->filename,ncatom);
    mprintf("         match those in associated parmtop (%i)!\n",P->natom);
    return 1;
  }
  stop=1;
  Frames=1; 
  close();
  return 0;
}

/*
 * AmberRestartNC::SetupWrite()
 * Setting up is done for each frame.
 */
int AmberRestartNC::SetupWrite() {
  return 0;
}

/*
 * AmberRestartNC::setupWriteForSet()
 * Create Netcdf restart file for the given frame and set it up. 
 */
int AmberRestartNC::setupWriteForSet(int set) {
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[2], count[2];
  char buffer[1024];
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };

  // Create filename for this set
  sprintf(buffer,"%s.%i",File->filename,set);

  // Create file
  if (checkNCerr(nc_create(buffer,NC_64BIT_OFFSET,&ncid),
    "Creating Netcdf restart file %s",buffer)) return 1;
  //if (debug>0) 
    mprintf("    Successfully created Netcdf restart file %s, ncid %i\n",buffer,ncid);

  // Time variable
  if (checkNCerr(nc_def_var(ncid,NCTIME,NC_DOUBLE,0,dimensionID,&timeVID),
    "Defining time variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,timeVID,"units",10,"picosecond"),
    "Writing time VID units.")) return 1;

  // Spatial dimension and variable
  if (checkNCerr(nc_def_dim(ncid,NCSPATIAL,3,&spatialDID),
    "Defining spatial dimension.")) return 1;
  dimensionID[0] = spatialDID;
  if (checkNCerr(nc_def_var(ncid,NCSPATIAL,NC_CHAR,1,dimensionID,&spatialVID),
    "Defining spatial variable.")) return 1;

  // Atom dimension
  if (checkNCerr(nc_def_dim(ncid,NCATOM,P->natom,&atomDID),
    "Defining atom dimension.")) return 1;
  ncatom=P->natom;

  // Coord variable
  dimensionID[0] = atomDID;
  dimensionID[1] = spatialDID;
  if (checkNCerr(nc_def_var(ncid,NCCOORDS,NC_DOUBLE,2,dimensionID,&coordVID),
    "Defining coordinates variable.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,coordVID,"units",8,"angstrom"),
    "Writing coordinates variable units.")) return 1;

  // Velocity variable
  if (hasVelocity) {
    if (checkNCerr(nc_def_var(ncid,NCVELO,NC_DOUBLE,2,dimensionID,&velocityVID),
      "Defining velocities variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,velocityVID,"units",19,"angstrom/picosecond"),
      "Writing velocities variable units.")) return 1;
    if (checkNCerr(nc_put_att_double(ncid,velocityVID,"scale_factor",NC_DOUBLE,1, &velocityScale),
      "Writing velocities scale factor.")) return 1;
  }

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
  if (BoxType!=0) {
    dimensionID[0]=cell_spatialDID;
    if (checkNCerr(nc_def_var(ncid,"cell_lengths",NC_DOUBLE,1,dimensionID,&cellLengthVID),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellLengthVID,"units",8,"angstrom"),
    "Writing cell length variable units.")) return 1;
    dimensionID[0]=cell_angularDID;
    if (checkNCerr(nc_def_var(ncid,"cell_angles",NC_DOUBLE,1,dimensionID,&cellAngleVID),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(ncid,cellAngleVID,"units",6,"degree"),
    "Writing cell angle variable units.")) return 1;
  }

  // Set up title
  if (title==NULL) {
    title=(char*) malloc(30*sizeof(char));
    strcpy(title,"Cpptraj Generated Restart");
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
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"Conventions",12,"AMBERRESTART"),
    "Writing conventions.")) return 1;
  if (checkNCerr(nc_put_att_text(ncid,NC_GLOBAL,"ConventionVersion",3,"1.0"),
    "Writing conventions version.")) return 1;

  // Replica temperature 
  if (hasTemperature) {
    mprintf("NETCDF: Defining replica temperature in output restart.\n");
    if ( checkNCerr(nc_def_var(ncid,NCTEMPERATURE,NC_DOUBLE,0,dimensionID,&TempVID),
         "Defining replica temperature in netcdf restart.") ) return 1;
    if ( checkNCerr(nc_put_att_text(ncid,TempVID,"units",6,"kelvin"),
         "Defining replica temperature units in netcdf restart.") ) return 1;
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

  // Close the file. It will be reopened write
  // NOTE: Do not close here. Since this is called for every frame a write immediately
  //       follows. Close in writeFrame.
  //close();
  
  return 0;
}


/*
 * AmberRestartNC::getFrame()
 * Get the specified frame from amber netcdf file
 * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
 */
int AmberRestartNC::getFrame(int set) {
  size_t start[2], count[2];

  // Get temperature
  if (TempVID!=-1) {
    if ( checkNCerr(nc_get_var_double(ncid, TempVID, &(F->T)),
                    "Getting replica temperature.")!=0 ) return 1;
    mprintf("DEBUG: Replica Temperature %lf\n",F->T);
  }

  // Read Coords 
  start[0]=0;
  start[1]=0;
  count[0]=ncatom;
  count[1]=3;
  if ( checkNCerr(nc_get_vara_double(ncid, coordVID, start, count, F->X),
                  "Getting Coords")!=0 ) return 1;

  // Read Velocity
  if (hasVelocity) {
    if (F->V==NULL) F->V = new Frame(ncatom,NULL);
    if ( checkNCerr(nc_get_vara_double(ncid, velocityVID, start, count, F->V->X),
                    "Getting velocities")!=0 ) return 1;
  }

  // Read box info 
  if (BoxType!=0) {
    count[0]=3;
    count[1]=0;
    if ( checkNCerr(nc_get_vara_double(ncid, cellLengthVID, start, count, F->box),
                    "Getting cell lengths.")!=0 ) return 1;
    if ( checkNCerr(nc_get_vara_double(ncid, cellAngleVID, start, count, &(F->box[3])),
                    "Getting cell angles.")!=0 ) return 1;
  }


  return 0;
}

/*
 * AmberRestartNC::writeFrame() 
 */
int AmberRestartNC::writeFrame(int set) {
  size_t start[2], count[2];

  // Set up file for this set
  if ( this->setupWriteForSet(set+OUTPUTFRAMESHIFT) ) return 1;

  // write coords
  start[0]=0;
  start[1]=0;
  count[0]=P->natom; // OK to use P->natom here since file was just set up
  count[1]=3;
  if (checkNCerr(nc_put_vara_double(ncid,coordVID,start,count,F->X),
      "Netcdf restart Writing coordinates %i",set)) return 1;

  // write velocity
  if (hasVelocity && F->V!=NULL) {
    if (checkNCerr(nc_put_vara_double(ncid,velocityVID,start,count,F->V->X),
        "Netcdf restart writing velocity %i",set)) return 1;
  }

  // write box
  if (BoxType!=0 && cellLengthVID!=-1) {
    count[0]=3;
    count[1]=0;
    if (checkNCerr(nc_put_vara_double(ncid,cellLengthVID,start,count,F->box),
      "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(ncid,cellAngleVID,start,count, &(F->box[3])),
      "Writing cell angles.")) return 1;
  }
  
  nc_sync(ncid); // Necessary after every write??

  // Close file for this set
  close();

  return 0;
}  

/*
 * Info()
 */
void AmberRestartNC::Info() {
  mprintf("is a NetCDF AMBER restart file");

  if (velocityVID!=-1) mprintf(", with velocities");
 
  if (TempVID!=-1) mprintf(", with replica temperature");

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
