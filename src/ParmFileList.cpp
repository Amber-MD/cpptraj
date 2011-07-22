// ParmList
#include <cstdlib>
#include <cstring>
#include "ParmFileList.h"
#include "CpptrajStdio.h"

// Constructor
ParmFileList::ParmFileList() {
  ParmList=NULL;
  Nparm=0;
  debug=0;
  hasCopies=false;
}

// Destructor
ParmFileList::~ParmFileList() {
  int i;

  if (ParmList!=NULL) {
    if (!hasCopies) {
      for (i=0; i<Nparm; i++) delete ParmList[i];
    }
    free(ParmList);
  }
}

// Set debug level
void ParmFileList::SetDebug(int debugIn) {
  if (debugIn>0) mprintf("ParmFileList debug level set to %i\n",debugIn);
  debug=debugIn;
}

// Return the parm structure with index num
AmberParm *ParmFileList::GetParm(int num) {
  if (num>=Nparm || num<0) return NULL;
  return ParmList[num];
}

/* ParmFileList::GetParm()
 * Return the parm structure based on arguments in the given arg list. 
 *   parm <parm name>
 *   parmindex <parm index>
 */
AmberParm *ParmFileList::GetParm(ArgList *A) {
  char *parmfilename;
  int pindex;
  AmberParm *P;
  // Get any parm keywords if present
  P=NULL;
  parmfilename=A->getKeyString("parm", NULL);
  pindex=A->getKeyInt("parmindex",0);
  // Associate trajectory with parameter file. Associate with default parm if none specified
  if (parmfilename!=NULL)
    pindex = this->GetParmIndex(parmfilename);
  P = this->GetParm(pindex);
  if (P==NULL) {
    mprinterr("    Error: Could not get parameter file:\n");
    mprinterr("           parmfilename=%s, pindex=%i\n",parmfilename,pindex);
    return NULL;
  }

  return P;
}

/*
 * ParmFileList::GetParmIndex()
 * Return the index in ParmList of the given Parm name. Use either the full
 * path or the base filename.
 */
int ParmFileList::GetParmIndex(char *name) {
  int i;
  int pindex;

  pindex=-1;
  for (i=0; i<Nparm; i++)
    if ( strcmp(name,ParmList[i]->parmfileName)==0 ||
         strcmp(name,ParmList[i]->parmName)==0 ) {
      pindex=i;
      break;
    }

  return pindex;
}

/* ParmFileList::Add()
 * Add a parameter file to the parm file list.
 */
int ParmFileList::Add(char *filename) {
  AmberParm *P;

  // Dont let a list that has copies add a new file
  if (hasCopies) {
    mprintf("    Warning: Attempting to add parm %s to a list that already\n",filename);
    mprintf("             has copies of parm files. This should not occur.\n");
    mprintf("             Skipping.\n");
    return 0;
  }

  // Check if this file has already been loaded
  if (GetParmIndex(filename)!=-1) {
    mprintf("    Warning: Parm %s already loaded, skipping.\n",filename);
    return 1;
  }

  P = new AmberParm();
  P->SetDebug(debug);

  if (P->OpenParm(filename)) {
    mprintf("Error: Could not open parm %s\n",filename);
    delete P;
    return 1;
  }

  if (debug>0) mprintf("    PARAMETER FILE %i: %s\n",Nparm,filename);
  // pindex is used for quick identification of the parm file
  P->pindex=Nparm;
  ParmList=(AmberParm**) realloc(ParmList,(Nparm+1) * sizeof(AmberParm*));
  ParmList[Nparm]=P;
  Nparm++;
  return 0;
}

/* ParmFileList::Add()
 * Add an existing AmberParm to parm file list. Currently used to keep track
 * of parm files corresponding to frames in the reference frame list.
 */
int ParmFileList::Add(AmberParm *ParmIn) {
  // Set the hasCopies flag so we know not to try and delete these parms
  hasCopies=true;
  //P->pindex=Nparm; // pindex should already be set
  ParmList=(AmberParm**) realloc(ParmList,(Nparm+1) * sizeof(AmberParm*));
  ParmList[Nparm]=ParmIn;
  Nparm++;
  return 0;
}

/*
 * ParmFileList::Replace()
 * Replace parm file at given position with newParm. If this list has only
 * copies do not delete the old parm, just replace.
 */
int ParmFileList::Replace(int num, AmberParm *newParm) {
  if (num>=Nparm || num<0) return 1;
  if (!hasCopies) delete ParmList[num];
  ParmList[num]=newParm;
  return 0;
}

// Print list of loaded parameter files
void ParmFileList::Print() {
  int i;

  mprintf("\nPARAMETER FILES:\n");
  if (Nparm==0) {
    mprintf("  No parameter files defined.\n");
    return;
  }

  for (i=0; i<Nparm; i++) {
    ParmList[i]->ParmInfo();
    //mprintf("  %i: %s, %i atoms (%i trajectory frames associated)\n",
    //        i,ParmList[i]->File.filename, ParmList[i]->natom, ParmList[i]->parmFrames);
  }
}
