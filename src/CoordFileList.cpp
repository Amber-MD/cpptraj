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
  trajfilename=NULL; 
  P = NULL;
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
 * CoordFileList::ProcessArgList()
 * Process arguments common to all trajectory lists. Currently gets traj 
 * filename and associated parm file.
 */
int CoordFileList::ProcessArgList(ArgList *A, ParmFileList *parmFileList) {
  char *parmfilename;
  int pindex;

  // Filename should be first arg
  trajfilename = A->getNextString();

  // Get any parm keywords if present
  parmfilename=A->getKeyString("parm", NULL);
  pindex=A->getKeyInt("parmindex",0);
  // Associate trajectory with parameter file. Associate with default parm if none specified
  if (parmfilename!=NULL)
    pindex=parmFileList->GetParmIndex(parmfilename);
  P = parmFileList->GetParm(pindex);
  if (P==NULL) {
    mprintf("    Error: Could not associate %s with a parameter file.\n",trajfilename);
    mprintf("    parmfilename=%s, pindex=%i\n",parmfilename,pindex);
    return 1;
  }
  if (debug>0)
    mprintf("    Associating traj %s with parm %s\n",trajfilename, P->parmName);

  return 0;
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

  basicTraj=new PtrajFile();
  // Set up basic file to determine (read) or set (write) type and format
  if (basicTraj->SetupFile(trajfilenameIN,fileAccess,writeFormat,writeType,debug)) {
    delete basicTraj;
    if (debug>1)
      mprintf("    Error: Could not set up file %s.\n",trajfilenameIN);
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
      mprintf("Error: SetupTrajectory(%s):\n",trajfilename);
      mprintf("       Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      T = new AmberRestartNC();
#else
      mprintf("Error: SetupTrajectory(%s):\n",trajfilename);
      mprintf("       Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case PDBFILE     : T = new PDBfile();      break;
    case CONFLIB     : T = new Conflib();      break;
    case MOL2FILE    : T = new Mol2File();     break;
    default:
      mprintf("    Error: Could not determine trajectory file %s type, skipping.\n",
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
  for (it = this->begin(); it != this->end(); it++)
    if ( strcmp(filenameIn, (*it)->File->filename)==0 ) return 1;

  return 0;
}

/*
 * CoordFileList::Info() - Call PrintInfo for each traj in the list.
 */
// May not be necessary
void CoordFileList::Info() {
  if (this->empty()) 
    mprintf("No files.\n");
  for (it = this->begin(); it != this->end(); it++)
    (*it)->PrintInfo(0);
}

