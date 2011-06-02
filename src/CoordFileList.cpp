// CoordFileList
#include <cstring> // strcmp
#include "CoordFileList.h"
#include "CpptrajStdio.h"
// All trajectory classes go here
#include "AmberTraj.h"
#ifdef BINTRAJ
#  include "AmberNetcdf.h"
#  include "AmberRestartNC.h"
#endif
#include "PDBfile.h"
#include "AmberRestart.h"
#include "RemdTraj.h"
#include "Conflib.h"
#include "Mol2File.h"

// CONSTRUCTOR
CoordFileList::CoordFileList() {
  debug=0;
  //trajfilename=NULL; 
  //P = NULL;
}

// DESTRUCTOR
CoordFileList::~CoordFileList() {
  //fprintf(stderr,"CoordFileList destructor\n");
  for (it = this->begin(); it != this->end(); it++)
    delete *it;
}

/*
 * CoordFileList::SetDebug()
 * Set trajectory list debug level.
 */
void CoordFileList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("CoordFileList() DEBUG LEVEL SET TO %i\n",debug);
}

/*
 * CoordFileList::SetupTrajectory()
 * Return a trajectory class with the File object set up for the specified
 * access; Associate with the given parm.
 */
TrajFile *CoordFileList::SetupTrajectory(char *trajfilenameIN, AccessType fileAccess,
                                        FileFormat writeFormat, FileType writeType) {
  PtrajFile *basicTraj;
  TrajFile *T;

  // Must specify a filename
  if (trajfilenameIN==NULL) {
    mprinterr("Error: CoordFileList::SetupTrajectory: No filename given.\n");
    return NULL;
  }

  basicTraj=new PtrajFile();
  // Set up basic file to determine (read) or set (write) type and format
  if (basicTraj->SetupFile(trajfilenameIN,fileAccess,writeFormat,writeType,debug)) {
    delete basicTraj;
    if (debug>1)
      mprinterr("    Error: Could not set up file %s.\n",trajfilenameIN);
    return NULL;
  }

  // Allocate specific trajectory type 
  switch ( basicTraj->fileFormat ) {
    case AMBERRESTART: T = new AmberRestart(); break;
    case AMBERTRAJ   : T = new AmberTraj();    break;
    case AMBERNETCDF : 
#ifdef BINTRAJ
      T = new AmberNetcdf();  
#else
      mprinterr("Error: SetupTrajectory(%s):\n",trajfilename);
      mprinterr("       Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      T = new AmberRestartNC();
#else
      mprinterr("Error: SetupTrajectory(%s):\n",trajfilename);
      mprinterr("       Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case PDBFILE     : T = new PDBfile();      break;
    case CONFLIB     : T = new Conflib();      break;
    case MOL2FILE    : T = new Mol2File();     break;
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",
              basicTraj->filename);
      delete basicTraj;
      return NULL;
  }
  // Set debug level
  T->debug=debug;
  // Place Ptraj File in Traj File
  T->File = basicTraj;
  // Set trajectory filename
  T->trajfilename = basicTraj->basefilename;
  //mprintf("DEBUG: trajfilename is %s\n",T->trajfilename);
  return T;
}


/*
 * CoordFileList::CheckFilename()
 * Check if filenameIn is already in use, return 1 if so.
 */
int CoordFileList::CheckFilename(char *filenameIn) {
  if (filenameIn==NULL) {
    mprinterr("Error: CoordFileList::CheckFilename: Called with NULL filename.\n");
    return 1;
  }
  for (it = this->begin(); it != this->end(); it++)
    if ( strcmp(filenameIn, (*it)->File->filename)==0 ) return 1;

  return 0;
}

/* CoordFileList::Info() - Call PrintInfo for each traj in the list.
 */
void CoordFileList::Info(int indent) {
  if (this->empty()) 
    mprintf("  No files.\n");
  for (it = this->begin(); it != this->end(); it++) {
    if (indent>0) mprintf("%*s",indent,"");
    (*it)->PrintInfo(0);
  }
}

