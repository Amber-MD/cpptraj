// ParmList
#include <cstring> // strcmp
#include "ParmFileList.h"
#include "CpptrajStdio.h"
#include "AtomMask.h"
#include "ParmFile.h"

// CONSTRUCTOR 
ParmFileList::ParmFileList() {
  Nparm=0;
  debug=0;
  hasCopies=false;
  bondsearch=false;
  molsearch=false;
}

// DESTRUCTOR
ParmFileList::~ParmFileList() {
    if (!hasCopies) {
      for (int i=0; i<Nparm; i++) delete ParmList[i];
    }
}

// ParmFileList::SetDebug()
/** Set debug level. */
void ParmFileList::SetDebug(int debugIn) {
  if (debugIn>0) mprintf("ParmFileList debug level set to %i\n",debugIn);
  debug=debugIn;
}

// ParmFileList::CheckCommand()
/** Check if the command in the arglist pertains to topology files.
  * \return 0 if command was recognized, 1 if not.
  */
int ParmFileList::CheckCommand(ArgList *argIn) {
  AtomMask tempMask;
  int pindex;
  // parm <filename> [<tag>]: Add <filename> to parm list
  if (argIn->CommandIs("parm")) {
    std::string parmtag = argIn->getNextTag();
    this->AddParmFile(argIn->getNextString(),parmtag);
    return 0;
  }
  // parmlist: Print list of loaded parm files
  if (argIn->CommandIs("parmlist")) {
    this->Print();
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
        ParmList[pindex]->SetupCharMask( tempMask, NULL );
        for (int atom=0; atom < ParmList[pindex]->natom; atom++) 
          if (tempMask.AtomInCharMask(atom)) ParmList[pindex]->AtomInfo(atom);
      } else {
        ParmList[pindex]->Summary();
      }
    } else
      mprinterr("Error: parminfo: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmwrite out <filename> [<parmindex>]: Write parm <parmindex> to <filename>
  if (argIn->CommandIs("parmwrite")) {
    char *outfilename = argIn->getKeyString("out",NULL);
    if (outfilename==NULL) {
      mprinterr("Error: parmwrite: No output filename specified (use 'out <filename>').\n");
      return 0;
    }
    pindex = argIn->getNextInteger(0);
    if (pindex < 0 || pindex >= Nparm) {
      mprinterr("Error: parmwrite: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    mprintf("\tWriting parm %i (%s) to Amber parm %s\n",pindex,
            ParmList[pindex]->parmName,outfilename);
    ParmFile pfile;
    pfile.SetDebug( debug );
    pfile.Write( *ParmList[pindex], outfilename, AMBERPARM );
    //ParmList[pindex]->WriteAmberParm(outfilename);
    return 0;
  }
  // parmstrip <mask> [<parmindex>]: Strip atoms int mask from parm
  if (argIn->CommandIs("parmstrip")) {
    char *mask0 = argIn->getNextMask();
    pindex = argIn->getNextInteger(0);
    if (pindex < 0 || pindex >= Nparm) {
      mprinterr("Error: parmstrip: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    tempMask.SetMaskString(mask0);
    // Since want to keep atoms outside mask, invert selection
    tempMask.InvertMask();
    ParmList[pindex]->SetupIntegerMask( tempMask, NULL );
    mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(), 
             ParmList[pindex]->natom - tempMask.Nselected, ParmList[pindex]->parmName);
    AmberParm *tempParm = ParmList[pindex]->modifyStateByMask(tempMask.Selected, NULL);
    if (tempParm==NULL) 
      mprinterr("Error: parmstrip: Could not strip parm.\n");
    else {
      tempParm->ParmInfo(ParmTags[pindex]);
      ReplaceParm(pindex, tempParm);
    }
    return 0;
  }
  // [parm]box [<parmindex>] [x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]
  //           [nobox]
  // Set the given parm box info to what is specified. If nobox, remove box info.
  if (argIn->CommandIs("box") || argIn->CommandIs("parmbox")) {
    double parmbox[6];
    parmbox[0] = argIn->getKeyDouble("x",0);
    parmbox[1] = argIn->getKeyDouble("y",0);
    parmbox[2] = argIn->getKeyDouble("z",0);
    parmbox[3] = argIn->getKeyDouble("alpha",0);
    parmbox[4] = argIn->getKeyDouble("beta",0);
    parmbox[5] = argIn->getKeyDouble("gamma",0);
    if (argIn->hasKey("nobox")) {
      parmbox[0]=-1; parmbox[1]=-1; parmbox[2]=-1;
      parmbox[3]=-1; parmbox[4]=-1; parmbox[5]=-1;
    }
    pindex = argIn->getNextInteger(0);
    if (pindex < 0 || pindex >= Nparm) {
      mprinterr("Error: box: parm index %i out of bounds.\n",pindex);
      return 0;
    }
    // Fill in missing parm box information from specified parm
    if (parmbox[0]==0) parmbox[0] = ParmList[pindex]->Box[0];
    if (parmbox[1]==0) parmbox[1] = ParmList[pindex]->Box[1];
    if (parmbox[2]==0) parmbox[2] = ParmList[pindex]->Box[2];
    if (parmbox[3]==0) parmbox[3] = ParmList[pindex]->Box[3];
    if (parmbox[4]==0) parmbox[4] = ParmList[pindex]->Box[4];
    if (parmbox[5]==0) parmbox[5] = ParmList[pindex]->Box[5];
    // Determine box type from parmbox angles
    ParmList[pindex]->boxType = CheckBoxType(parmbox+3, 1);
    ParmList[pindex]->Box[0] = parmbox[0];
    ParmList[pindex]->Box[1] = parmbox[1];
    ParmList[pindex]->Box[2] = parmbox[2];
    ParmList[pindex]->Box[3] = parmbox[3];
    ParmList[pindex]->Box[4] = parmbox[4];
    ParmList[pindex]->Box[5] = parmbox[5];
    return 0;
  } 
  // parmbondinfo [<parmindex>]: Print bond information for parm <parmindex>
  //     (0 by default).
  if (argIn->CommandIs("parmbondinfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm) 
      ParmList[pindex]->PrintBondInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmmolinfo [<parmindex>]: Print molecule information for parm
  //     <parmindex> (0 by default).
  if (argIn->CommandIs("parmmolinfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm)
      ParmList[pindex]->PrintMoleculeInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // parmresinfo [<parmindex>]: Print residue information for parm
  //     <parmindex> (0 by default).
  if (argIn->CommandIs("parmresinfo")) {
    pindex = argIn->getNextInteger(0);
    if (pindex>=0 && pindex<Nparm)
      ParmList[pindex]->PrintResidueInfo();
    else
      mprinterr("Error: parm %i not loaded.\n",pindex);
    return 0;
  }
  // bondsearch: Indicate that if bond information not found in topology
  //     it should be determined by distance search.
  if (argIn->CommandIs("bondsearch")) {
    mprintf("\tInfo: Bond info will be determined from distance search if not present.\n");
    bondsearch=true;
    return 0;
  }
  // molsearch: Indicate that if molecule information not found in 
  //     topology file it should be determined by bonding information.
  if (argIn->CommandIs("molsearch")) {
    mprintf("\tInfo: Molecule info will be determined from bonds if not present.\n");
    molsearch=true;
    return 0;
  }
  // nobondsearch: Turn off bond search.
  if (argIn->CommandIs("nobondsearch")) {
    mprintf("\tInfo: Bond search is off.\n");
    bondsearch=false;
    return 0;
  }
  // nomolsearch: Turn off molecule search.
  if (argIn->CommandIs("nomolsearch")) {
    mprintf("\tInfo: Molecule search is off.\n");
    molsearch=false;
    return 0;
  }
  // Unrecognized parm command
  return 1;
}

// ParmFileList::GetParm()
/** Return the parm structure with index num. */
AmberParm *ParmFileList::GetParm(int num) {
  if (num>=Nparm || num<0) return NULL;
  return ParmList[num];
}

// ParmFileList::GetParm()
/** Return the parm structure based on arguments in the given arg list. 
  *   parm <parm name>
  *   parmindex <parm index>
  * \param argIn argument list that contains parm-related keyword
  * \return parm specified by 'parm' or 'parmindex', or the first parm. NULL on error.
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

// ParmFileList::GetParmIndex()
/** Return the index in ParmList of the given Parm name. Use either the full
  * path, the base filename, or a tag.
  */
int ParmFileList::GetParmIndex(char *name) {
  int pindex;
  std::string ParmTag;

  // if first char of name is a bracket, assume tag.
  if (name[0]=='[') {
    ParmTag.assign( name );
    pindex = GetParmIndexByTag( ParmTag );
  
  // Otherwise assume filename or base filename
  } else {
    pindex=-1;
    for (int i=0; i<Nparm; i++)
      if ( strcmp(name,ParmList[i]->parmfileName)==0 ||
           strcmp(name,ParmList[i]->parmName)==0 ) {
        pindex=i;
        break;
      }
  }

  return pindex;
}

// ParmFileList::GetParmIndexByTag()
/** Return index of parm that matches tag. */
int ParmFileList::GetParmIndexByTag(std::string &ParmTag) {
  int i; 
  if (ParmTag.empty()) return -1;
  for (i = 0; i < Nparm; i++)
    if ( ParmTags[i].compare( ParmTag )==0 ) return i;
  return -1;
}

// ParmFileList::AddParmFile()
/** Add a parameter file to the parm file list. */
int ParmFileList::AddParmFile(char *filename) {
  std::string emptystring;
  return AddParmFile(filename, emptystring);
}

// ParmFileList::AddParmFile()
/** Add a parameter file to the parm file list with optional tag. */
int ParmFileList::AddParmFile(char *filename, std::string &ParmTag) {
  AmberParm *P;
  ParmFile pfile;

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

  // If tag specified, check if tag already in use
  if (GetParmIndexByTag(ParmTag)!=-1) {
    mprintf("    Warning: Parm tag [%s] already in use.\n",ParmTag.c_str());
    return 1;
  }

  P = new AmberParm();
  P->SetDebug(debug);
  pfile.SetDebug( debug );
  int err = pfile.Read(*P, filename, bondsearch, molsearch);
  if (err!=0) {
  //if (P->OpenParm(filename,bondsearch,molsearch)) {
    mprinterr("Error: Could not open parm %s\n",filename);
    delete P;
    return 1;
  }

  if (debug>0) mprintf("    PARAMETER FILE %i: %s\n",Nparm,filename);
  // pindex is used for quick identification of the parm file
  P->pindex=Nparm;
  ParmList.push_back(P);
  ParmTags.push_back( ParmTag );
  ++Nparm;
  return 0;
}

// ParmFileList::AddParm()
/** Add an existing AmberParm to parm file list. Currently used to keep track
  * of parm files corresponding to frames in the reference frame list.
  */
int ParmFileList::AddParm(AmberParm *ParmIn) {
  // Set the hasCopies flag so we know not to try and delete these parms
  hasCopies=true;
  //P->pindex=Nparm; // pindex should already be set
  ParmList.push_back(ParmIn);
  Nparm++;
  return 0;
}

// ParmFileList::ReplaceParm()
/** Replace parm file at given position with newParm. If this list has only
  * copies do not delete the old parm, just replace.
  */
int ParmFileList::ReplaceParm(int num, AmberParm *newParm) {
  if (num>=Nparm || num<0) return 1;
  if (!hasCopies) delete ParmList[num];
  ParmList[num]=newParm;
  return 0;
}

// ParmFileList::Print()
/// Print list of loaded parameter files
void ParmFileList::Print() {
  mprintf("\nPARAMETER FILES:\n");
  if (Nparm==0) {
    mprintf("  No parameter files defined.\n");
    return;
  }

  for (int i=0; i<Nparm; i++) {
    ParmList[i]->ParmInfo(ParmTags[i]);
    //mprintf("  %i: %s, %i atoms (%i trajectory frames associated)\n",
    //        i,ParmList[i]->File.filename, ParmList[i]->natom, ParmList[i]->parmFrames);
  }
}
