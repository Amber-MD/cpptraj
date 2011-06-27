// PDBfile
#include "Traj_PDBfile.h"
#include <cstdlib>
#include "PDBfileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
PDBfile::PDBfile() {
  pdbAtom=0;
  pdbWriteMode=SINGLE;
  dumpq=false;
  pdbAtomNames=NULL;

  trajResNames=NULL;
  trajAtomsPerMol=NULL;
  trajResNums=NULL;
  trajCharges=NULL;
  trajRadii=NULL;
}

// DESTRUCTOR
PDBfile::~PDBfile() { 
  //fprintf(stderr,"PDBfile Destructor.\n");
  // If WRITE/APPEND pdbAtomNames were passed in from a parm so dont free
  if (tfile->access==READ && pdbAtomNames!=NULL) 
    free(pdbAtomNames);
}
//------------------------------------------------------------------------
/* PDBfile::closeTraj()
 */
void PDBfile::closeTraj() {
  // On WRITE only close if not writing 1 pdb per frame
  if (tfile->access==WRITE) {
    if ( pdbWriteMode!=MULTI) {
      // Dont attempt a write if not successfully set up
      //if (skip==1)
      if (tfile->IsOpen()) tfile->IO->Printf("%-6s\n","END");
      tfile->CloseFile();
    }
  // Otherwise just close the file
  } else
    tfile->CloseFile();
}

/* PDBfile::openTraj()
 */
int PDBfile::openTraj() {
  int err;

  err = 0; 
  switch (tfile->access) {
    case READ  : err = tfile->OpenFile(); break;
    case WRITE :
      // If writing 1 pdb per frame do not open here
      if (pdbWriteMode!=MULTI) err = tfile->OpenFile(); 
      break;
    case APPEND:
      // NOTE: Test this, should actually be ok
      mprinterr("Error: %s: Append currently not supported for PDB files.\n",tfile->filename);
      err=1;
      break;
  }
  
  return err;
}

/* PDBfile::setupRead()
 * Scan PDB file to determine number of frames (models). The first frame will 
 * also be checked to ensure that the atom names match those in the parm file
 * in TrajectoryFile.
 */
int PDBfile::setupRead(int natom) {
  int atom, Frames;
  bool scanPDB;

  if ( this->openTraj() ) return -1;

  // Allocate space for checking pdb atom names
  pdbAtomNames = (NAME*) realloc(pdbAtomNames,natom * sizeof(NAME));
  // Two strats - check for MODEL keywords or see how many times natom ATOMs can be read
  // Currently employing the latter.
  Frames=0;
  scanPDB=true;
  while (scanPDB) {
    atom=0;
    while ( atom < natom ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( tfile->IO->Gets(buffer,256) ) { scanPDB=false; break; }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (!isPDBatomKeyword(buffer)) continue;
      // If still on frame 1, store atom names in order to subsequently check
      // if they match those in associated parm (TrajectoryFile)
      if (Frames==0) {
        pdb_name(buffer,(char*)(pdbAtomNames + atom));
      }
      atom++;
    }
    if (Frames==0) { 
      pdbAtom = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom!=atom) {
        mprinterr("Error: %s: Reading frame %i, got %i atoms, expected %i.\n",tfile->filename,
                  Frames,atom,pdbAtom);
        return -1;
      }
    }  
    if (scanPDB) Frames++;
  }
  this->closeTraj();

  if (Frames<1) {
    mprinterr("Error: %s: No frames read. atom=%i expected %i.\n",tfile->filename,
            atom,natom);
    return -1;
  }
  if (debug>0) mprintf("PDBfile: %s has %i atoms, %i frames.\n",tfile->filename,
                       pdbAtom,Frames);
  return Frames;
}

/* PDBfile::CheckPdbNames
 * Check previously read in PDB names (setupRead) against the provided
 * names, print a warning if they dont match up.
 * Return the number of mismatches
 */
int PDBfile::CheckPdbNames(NAME *names) {
  int numMismatch = 0;
  if (pdbAtomNames==NULL) return -1;
  for (int atom=0; atom < pdbAtom; atom++) {
    if ( (names[atom][0] != pdbAtomNames[atom][0]) ||
         (names[atom][1] != pdbAtomNames[atom][1]) ||
         (names[atom][2] != pdbAtomNames[atom][2]) ||
         (names[atom][3] != pdbAtomNames[atom][3])   ) 
    { 
      mprintf("Warning: %s: Atom %i name [%s] does not match parm name [%s]\n",
              tfile->filename,atom,pdbAtomNames[atom],names[atom]);
      numMismatch++;
    }
  }
  return numMismatch;
}

/* PDBfile::readFrame()
 * Read frame (model) from PDB file. 
 */
int PDBfile::readFrame(int set,double *X, double *box, double *T) {
  int atom, atom3;

  atom=0;
  atom3=0;
  while (atom < pdbAtom) {
    if ( tfile->IO->Gets(buffer,256) ) return 1;
    // Skip non-ATOM records
    if (!isPDBatomKeyword(buffer)) continue;
    // Read current PDB record XYZ into Frame
    pdb_xyz(buffer,X+atom3);
    atom++; 
    atom3+=3;
  }

  return 0;
}

/* PDBfile::setupWrite()
 */ 
int PDBfile::setupWrite(int natom ) {
  pdbAtom = natom;
  // Check that SetParmInfo has indeed been called
  if (pdbAtomNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Atom names are NULL.\n",tfile->filename); 
    return 1;
  }
  if (trajResNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue names are NULL.\n",tfile->filename); 
    return 1;
  }
  if (trajAtomsPerMol==NULL) {
    mprinterr("Error: setupWrite [%s]: Atoms per molecule is NULL.\n",tfile->filename); 
    return 1;
  }
  if (trajResNums==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue #s are NULL.\n",tfile->filename); 
    return 1;
  }
  if (dumpq && trajCharges==NULL) {
    mprinterr("Error: setupWrite [%s]: Charges are NULL.\n",tfile->filename); 
    return 1;
  }
  return 0;
}

/* PDBfile::SetWriteMode()
 * Set write mode to SINGLE, MODEL, or MULTI
 */
void PDBfile::SetWriteMode(PDBWRITEMODE modeIn) { 
  pdbWriteMode = modeIn; 
  //mprintf("PDB WRITE MODE SET TO %i\n",(int)pdbWriteMode);
}

/* PDBfile::NumFramesToWrite()
 * Modify write mode based on expected number of frames being written.
 */
void PDBfile::NumFramesToWrite(int parmFrames) {
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use MODEL keyword to separate frames.
  if (pdbWriteMode==SINGLE && parmFrames>1) pdbWriteMode=MODEL;
}

/* PDBfile::SetParmInfo()
 * In order to write pdb files, need some parm info like atom names etc.
 * This should be called before setupWrite.
 */
void PDBfile::SetParmInfo(NAME *names, NAME *resnames, int *atomsPerMol,
                          int *resnums, double *charge, double *radii) {
  //mprintf("SETTING PARM INFO FOR PDB FILE\n"); // DEBUG
  pdbAtomNames = names;
  trajResNames = resnames;
  trajAtomsPerMol = atomsPerMol;
  trajResNums = resnums;
  trajCharges = charge;
  trajRadii = radii;
  // DEBUG
  //mprintf("Atomname 0 %s\n",trajAtomNames[0]);
  //mprintf("Resname 0  %s\n",trajResNames[0]);
  //mprintf("AtPerMol 0 %i\n",trajAtomsPerMol[0]);
  //mprintf("ResNums 0  %i\n",trajResNums[0]);
  //mprintf("Charge 0   %lf\n",trajCharges[0]);
}

/* PDBfile::writeFrame()
 * Write the frame (model) to PDB file.
 */
int PDBfile::writeFrame(int set,double *X,double *box,double T) {
  int i,i3,res,atom,mol,lastAtomInMol,bufferSize;
  float Occ, B;

  // If writing 1 pdb per frame set up output filename and open
  if (pdbWriteMode==MULTI) {
    NumberFilename(buffer,tfile->filename,set + OUTPUTFRAMESHIFT);
    if (tfile->IO->Open(buffer,"wb")) return 1;
  // If specified, write MODEL keyword
  } else if (pdbWriteMode==MODEL) {
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    tfile->IO->Printf("MODEL     %i\n",set + OUTPUTFRAMESHIFT);
  }

  res=0; Occ=0.0; B=0.0;
  i3=0;
  atom=1; // Actual PDB atom number
  mol=0;
  if (trajAtomsPerMol!=NULL)
    lastAtomInMol=trajAtomsPerMol[0];
  for (i=0; i<pdbAtom; i++) {
    // If this atom belongs to a new molecule print a TER card
    /*if (P->atomsPerMol!=NULL) {
      if (i == lastAtomInMol) {
        pdb_write_ATOM(buffer,PDBTER,atom,(char*)"",P->ResidueName(res),'X',res+1,
                       0,0,0,0,0,(char*)"\0",dumpq);
        tfile->IO->Write(buffer,sizeof(char),strlen(buffer));
        atom++;
        mol++;
        if (mol<P->molecules) lastAtomInMol += P->atomsPerMol[mol];
      }
    }*/
    // figure out the residue number
    if ( i==trajResNums[res+1] ) res++;
    if (dumpq) Occ = (float) trajCharges[i]; 
    bufferSize=pdb_write_ATOM(buffer,PDBATOM,atom,pdbAtomNames[i],trajResNames[res],'X',
                              res+1,X[i3],X[i3+1],X[i3+2],Occ,B,(char*)"\0",dumpq);
    tfile->IO->Write(buffer,sizeof(char),bufferSize); 
    i3+=3;
    atom++;
  }

  // If writing 1 pdb per frame, close output file
  if (pdbWriteMode==MULTI) {
    tfile->IO->Printf("%-6s\n","END");
    tfile->IO->Close();
  // If MODEL keyword was written, write corresponding ENDMDL record
  } else if (pdbWriteMode==MODEL) {
    tfile->IO->Printf("ENDMDL\n");
  }

  return 0;
}

/* PDBfile::info()
 */
void PDBfile::info() {
  mprintf("is a PDB file");
  if (tfile->access==WRITE) {
    if (pdbWriteMode==MULTI)
      mprintf(" (1 file per frame)");
    else if (pdbWriteMode==MODEL)
      mprintf(" (1 MODEL per frame)");
    if (dumpq) mprintf(", writing out charges to occupancy column");
  }
   
/*    if (p->option2 == 1) 
      printfone(" with no atom wrapping");
    if (p->option1 == 1)
      printfone(": AMBER charges and radii in prmtop to occupancy and temp factor columns");
    else if (p->option1 == 2)
      printfone(": AMBER charges and PARSE radii to occupancy and temp factor columns");
    else if (p->option1 == 3)
      printfone(": AMBER charges and vdw radii (r*) to occupancy and temp factor columns");*/
}
