// ParmList
#include <cstdlib>
#include <cstring>
#include "ParmFileList.h"
#include "CpptrajStdio.h"
#include "AtomMask.h"

// CONSTRUCTOR 
ParmFileList::ParmFileList() {
  ParmList=NULL;
  Nparm=0;
  debug=0;
  hasCopies=false;
  bondsearch=false;
  molsearch=false;
}

// DESTRUCTOR
ParmFileList::~ParmFileList() {
  if (ParmList!=NULL) {
    if (!hasCopies) {
      for (int i=0; i<Nparm; i++) delete ParmList[i];
    }
    free(ParmList);
  }
}

/* ParmFileList::SetDebug()
 * Set debug level.
 */
void ParmFileList::SetDebug(int debugIn) {
  if (debugIn>0) mprintf("ParmFileList debug level set to %i\n",debugIn);
  debug=debugIn;
}

/* ParmFileList::CheckCommand()
 * Check if the command in the arglist pertains to topology files.
 * Return 0 if command was recognized, 1 if not.
 */
int ParmFileList::CheckCommand(ArgList *argIn) {
  AtomMask tempMask;
  int pindex;
  // parm <filename>: Add <filename> to parm list
  if (argIn->CommandIs("parm")) {
    this->Add(argIn->getNextString());
    return 0;
  }
  // parminfo [<parmindex>] [<mask>]: Print information on parm <parmindex> 
  //     (0 by default). If <mask> is given print info on atoms in mask. If
  //     no mask given print overall information.
  if (argIn->CommandIs("parminfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm) {
      char *maskarg = argIn->getNextMask();
      if (maskarg!=NULL) {
        tempMask.SetMaskString( maskarg );
        tempMask.SetupCharMask( ParmList[pindex], debug);
        for (int atom=0; atom < ParmList[pindex]->natom; atom++) 
          if (tempMask.AtomInCharMask(atom)) ParmList[pindex]->AtomInfo(atom);
      } else {
        ParmList[pindex]->Summary();
      }
    } else
      mprinterr("\tError: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmbondinfo [<parmindex>]: Print bond information for parm <parmindex>
  //     (0 by default).
  if (argIn->CommandIs("parmbondinfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm) 
      ParmList[pindex]->PrintBondInfo();
    else
      mprinterr("\tError: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmmolinfo [<parmindex>]: Print molecule information for parm
  //     <parmindex> (0 by default).
  if (argIn->CommandIs("parmmolinfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm)
      ParmList[pindex]->PrintMoleculeInfo();
    else
      mprinterr("\tError: parm %i not loaded.\n",pindex);
    return 0;
  }
  // bondsearch: Indicate that if bond information not found in topology
  //     it should be determined by distance search.
  if (argIn->CommandIs("bondsearch")) {
    bondsearch=true;
    return 0;
  }
  // molsearch: Indicate that if molecule information not found in 
  //     topology file it should be determined by bonding information.
  if (argIn->CommandIs("molsearch")) {
    molsearch=true;
    return 0;
  }
  // nobondsearch: Turn off bond search.
  if (argIn->CommandIs("nobondsearch")) {
    bondsearch=false;
    return 0;
  }
  // nomolsearch: Turn off molecule search.
  if (argIn->CommandIs("nomolsearch")) {
    molsearch=false;
    return 0;
  }
  // Unrecognized parm command
  return 1;
}

/* ParmFileList::GetParm()
 * Return the parm structure with index num.
 */
AmberParm *ParmFileList::GetParm(int num) {
  if (num>=Nparm || num<0) return NULL;
  return ParmList[num];
}

/* ParmFileList::GetParm()
 * Return the parm structure based on arguments in the given arg list. 
 *   parm <parm name>
 *   parmindex <parm index>
 */
AmberParm *ParmFileList::GetParm(ArgList &argIn) {
  char *parmfilename;
  int pindex;
  AmberParm *P;
  // Get any parm keywords if present
  P=NULL;
  parmfilename=argIn.getKeyString("parm", NULL);
  pindex=argIn.getKeyInt("parmindex",0);
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

/* ParmFileList::GetParmIndex()
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

  if (P->OpenParm(filename,bondsearch,molsearch)) {
    mprinterr("Error: Could not open parm %s\n",filename);
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

/* ParmFileList::Replace()
 * Replace parm file at given position with newParm. If this list has only
 * copies do not delete the old parm, just replace.
 */
int ParmFileList::Replace(int num, AmberParm *newParm) {
  if (num>=Nparm || num<0) return 1;
  if (!hasCopies) delete ParmList[num];
  ParmList[num]=newParm;
  return 0;
}

/* ParmFileList::Print()
 * Print list of loaded parameter files
 */
void ParmFileList::Print() {
  mprintf("\nPARAMETER FILES:\n");
  if (Nparm==0) {
    mprintf("  No parameter files defined.\n");
    return;
  }

  for (int i=0; i<Nparm; i++) {
    ParmList[i]->ParmInfo();
    //mprintf("  %i: %s, %i atoms (%i trajectory frames associated)\n",
    //        i,ParmList[i]->File.filename, ParmList[i]->natom, ParmList[i]->parmFrames);
  }
}
