// PDBfile
#include <cstdlib>
#include "Traj_PDBfile.h"
#include "PDBfileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
PDBfile::PDBfile() {
  pdbAtom=0;
  pdbWriteMode=SINGLE;
  dumpq=false;
  dumpr=false;

  pdbAtomNames=NULL;
  trajResNames=NULL;
  trajAtomsPerMol=NULL;
  trajResNums=NULL;
  trajCharges=NULL;
  trajRadii=NULL;

  chainID=NULL;
  chainchar=' ';
}

// DESTRUCTOR
PDBfile::~PDBfile() { 
  //fprintf(stderr,"PDBfile Destructor.\n");
  if (chainID!=NULL) delete[] chainID;
}
//------------------------------------------------------------------------
// PDBfile::closeTraj()
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

// PDBfile::openTraj()
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

// PDBfile::setupRead()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int PDBfile::setupRead(AmberParm *trajParm) {
  int atom, Frames;
  int numMismatch = 0;
  bool scanPDB;
  char pdbAtomName[5];

  if ( this->openTraj() ) return -1;

  // Two strats - check for MODEL keywords or see how many times natom ATOMs can be read
  // Currently employing the latter.
  Frames=0;
  scanPDB=true;
  while (scanPDB) {
    atom=0;
    while ( atom < trajParm->natom ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( tfile->IO->Gets(buffer,256) ) { scanPDB=false; break; }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (!isPDBatomKeyword(buffer)) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        pdb_name(buffer,pdbAtomName);
        if (!trajParm->AtomNameIs(atom, pdbAtomName)) {
          if (debug>1) 
            mprintf("Warning: %s: Atom %i name [%s] does not match parm atom name [%s]\n",
                    tfile->filename,atom,pdbAtomName,trajParm->AtomName(atom));
          numMismatch++;
        }
      }
      atom++;
    }
    if (Frames==0) { 
      pdbAtom = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom!=atom) {
        mprintf("Warning: %s: Reading frame %i, got %i atoms, expected %i.\n",tfile->filename,
                  Frames+OUTPUTFRAMESHIFT,atom,pdbAtom);
        mprintf("         Only using frames 1-%i\n",Frames);
        scanPDB=false;
        break;
      }
    }  
    if (scanPDB) Frames++;
  }
  this->closeTraj();

  if (Frames<1) {
    mprinterr("Error: %s: No frames read. atom=%i expected %i.\n",tfile->filename,
            atom,trajParm->natom);
    return -1;
  }
  if (debug>0) mprintf("PDBfile: %s has %i atoms, %i frames.\n",tfile->filename,
                       pdbAtom,Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",tfile->filename,
            numMismatch,trajParm->parmName);

  return Frames;
}

// PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int PDBfile::readFrame(int set,double *X, double *V,double *box, double *T) {
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

// PDBfile::processWriteArgs()
int PDBfile::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("dumpq")) this->SetDumpq();
  if (argIn->hasKey("model")) this->SetWriteMode(MODEL);
  if (argIn->hasKey("multi")) this->SetWriteMode(MULTI);
  char *temp = argIn->getKeyString("chainid",NULL);
  if (temp!=NULL) chainchar = temp[0];
  return 0;
}

// PDBfile::setupWrite()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int PDBfile::setupWrite(AmberParm *trajParm) {
  pdbAtom = trajParm->natom;
  pdbAtomNames = trajParm->AtomNames_ptr();
  trajResNames = trajParm->ResidueNames_ptr();
  trajAtomsPerMol = trajParm->AtomsPerMol_ptr();
  trajResNums = trajParm->ResAtomNums_ptr();
  trajCharges = trajParm->Charges_ptr();
  trajRadii = trajParm->GB_radii_ptr();

  // Set a chainID for each atom
  if (chainID!=NULL) delete[] chainID;
  chainID = new char[pdbAtom];
  for (int atom=0; atom < pdbAtom; atom++) chainID[atom] = chainchar;

  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (pdbWriteMode==SINGLE && trajParm->parmFrames>1) pdbWriteMode=MODEL;
  // Check that all parm info needed is present
  if (pdbAtomNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Atom names are NULL.\n",tfile->filename); 
    return 1;
  }
  if (trajResNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue names are NULL.\n",tfile->filename); 
    return 1;
  }
  if (trajAtomsPerMol==NULL) {
    if (debug>0) {
      mprintf("Warning: setupWrite [%s]: Atoms per molecule is NULL.\n",tfile->filename); 
      mprintf("         TER cards will not be written to PDB.\n");
    }
  }
  if (trajResNums==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue #s are NULL.\n",tfile->filename); 
    return 1;
  }
  if (dumpq && trajCharges==NULL) {
    mprinterr("Error: setupWrite [%s]: Charges are NULL.\n",tfile->filename); 
    return 1;
  }
  if (dumpr && trajRadii==NULL) {
    mprintf("Warning: setupWrite[%s]: Radii are NULL and will not be printed.\n",
            tfile->filename);
    dumpr=false;
  }
  return 0;
}

// PDBfile::SetWriteMode()
/** Set write mode to SINGLE, MODEL, or MULTI */
void PDBfile::SetWriteMode(PDBWRITEMODE modeIn) { 
  pdbWriteMode = modeIn; 
  //mprintf("PDB WRITE MODE SET TO %i\n",(int)pdbWriteMode);
}

// PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int PDBfile::writeFrame(int set,double *X,double *V,double *box,double T) {
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
  lastAtomInMol=-1;
  if (trajAtomsPerMol!=NULL)
    lastAtomInMol=trajAtomsPerMol[0];
  for (i=0; i<pdbAtom; i++) {
    // If this atom belongs to a new molecule print a TER card
    if (i == lastAtomInMol) {
      bufferSize=pdb_write_ATOM(buffer,PDBTER,atom,(char*)"",trajResNames[res],chainID[i],res+1,
                                0,0,0,0,0,(char*)"\0",dumpq);
      tfile->IO->Write(buffer,sizeof(char),bufferSize);
      atom++;
      mol++;
      lastAtomInMol += trajAtomsPerMol[mol];
    }
    // figure out the residue number
    if ( i==trajResNums[res+1] ) res++;
    if (dumpq) Occ = (float) trajCharges[i];
    if (dumpr) B = (float) trajRadii[i]; 
    bufferSize=pdb_write_ATOM(buffer,PDBATOM,atom,pdbAtomNames[i],trajResNames[res],chainID[i],
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

// PDBfile::info()
void PDBfile::info() {
  mprintf("is a PDB file");
  if (tfile->access==WRITE) {
    if (pdbWriteMode==MULTI)
      mprintf(" (1 file per frame)");
    else if (pdbWriteMode==MODEL)
      mprintf(" (1 MODEL per frame)");
    if (dumpq && !dumpr) 
      mprintf(", writing charges to occupancy column");
    else if (dumpr && !dumpq) 
      mprintf(", writing GB radii to B-factor column");
    else if (dumpr && dumpq)
      mprintf(", writing charges/GB radii to occupancy/B-factor columns");
  }
}
