// PDBfile
#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"

const size_t PDBfile::BUF_SIZE = 83;

// CONSTRUCTOR
PDBfile::PDBfile() {
  pdbAtom_=0;
  pdbWriteMode_=SINGLE;
  dumpq_=false;
  dumpr_=false;

  pdbAtomNames_=NULL;
  trajResNames_=NULL;
  trajAtomsPerMol_=NULL;
  trajResNums_=NULL;
  trajCharges_=NULL;
  trajRadii_=NULL;

  chainchar_=' ';
}

//------------------------------------------------------------------------
// PDBfile::closeTraj()
void PDBfile::closeTraj() {
  // On WRITE only close if not writing 1 pdb per frame
  if (access_==WRITE) {
    if ( pdbWriteMode_!=MULTI) {
      // Dont attempt a write if not successfully set up
      //if (skip==1)
      if (isOpen_) Printf("%-6s\n","END");
      CloseFile();
    }
  // Otherwise just close the file
  } else
    CloseFile();
}

// PDBfile::openTraj()
int PDBfile::openTraj() {
  int err;

  err = 0; 
  switch (access_) {
    case READ  : err = OpenFile(); break;
    case WRITE :
      // If writing 1 pdb per frame do not open here
      if (pdbWriteMode_!=MULTI) err = OpenFile(); 
      break;
    case APPEND:
      // NOTE: Test this, should actually be ok
      mprinterr("Error: %s: Append currently not supported for PDB files.\n",BaseName());
      err=1;
      break;
  }
  
  return err;
}

// PDBfile::setupTrajin()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int PDBfile::setupTrajin(Topology *trajParm) {
  char buffer[BUF_SIZE];
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
      if ( IO->Gets(buffer,BUF_SIZE) ) { scanPDB=false; break; }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (!isPDBatomKeyword(buffer)) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        pdb_name(buffer,pdbAtomName);
        if (!trajParm->AtomNameIs(atom, pdbAtomName)) {
          if (debug_>1) 
            mprintf("Warning: %s: Atom %i name [%s] does not match parm atom name [%s]\n",
                    BaseName(),atom,pdbAtomName,trajParm->AtomName(atom));
          numMismatch++;
        }
      }
      atom++;
    }
    if (Frames==0) { 
      pdbAtom_ = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom_!=atom) {
        mprintf("Warning: %s: Reading frame %i, got %i atoms, expected %i.\n",BaseName(),
                  Frames+OUTPUTFRAMESHIFT,atom,pdbAtom_);
        mprintf("         Only using frames 1-%i\n",Frames);
        scanPDB=false;
        break;
      }
    }  
    if (scanPDB) Frames++;
  }
  this->closeTraj();

  if (Frames<1) {
    mprinterr("Error: %s: No frames read. atom=%i expected %i.\n",BaseName(),
            atom,trajParm->natom);
    return -1;
  }
  if (debug_>0) mprintf("PDBfile: %s has %i atoms, %i frames.\n",BaseName(),
                       pdbAtom_,Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",BaseName(),
            numMismatch,trajParm->parmName);

  return Frames;
}

// PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int PDBfile::readFrame(int set,double *X, double *V,double *box, double *T) {
  char buffer[BUF_SIZE];
  int atom, atom3;

  atom=0;
  atom3=0;
  while (atom < pdbAtom_) {
    if ( IO->Gets(buffer,BUF_SIZE) ) return 1;
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
  if (temp!=NULL) chainchar_ = temp[0];
  return 0;
}

// PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int PDBfile::setupTrajout(Topology *trajParm) {
  pdbAtom_ = trajParm->natom;
  pdbAtomNames_ = trajParm->AtomNames_ptr();
  trajResNames_ = trajParm->ResidueNames_ptr();
  trajAtomsPerMol_ = trajParm->AtomsPerMol_ptr();
  trajResNums_ = trajParm->ResAtomNums_ptr();
  trajCharges_ = trajParm->Charges_ptr();
  trajRadii_ = trajParm->GB_radii_ptr();

  // Set a chainID for each atom
  chainID_.clear();
  chainID_.resize(pdbAtom_, chainchar_);

  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (pdbWriteMode_==SINGLE && trajParm->parmFrames>1) 
    pdbWriteMode_=MODEL;
  // Check that all parm info needed is present
  if (pdbAtomNames_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Atom names are NULL.\n",BaseName()); 
    return 1;
  }
  if (trajResNames_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Residue names are NULL.\n",BaseName()); 
    return 1;
  }
  if (trajAtomsPerMol_==NULL) {
    if (debug_>0) {
      mprintf("Warning: setupTrajout [%s]: Atoms per molecule is NULL.\n",BaseName()); 
      mprintf("         TER cards will not be written to PDB.\n");
    }
  }
  if (trajResNums_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Residue #s are NULL.\n",BaseName()); 
    return 1;
  }
  if (dumpq_ && trajCharges_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Charges are NULL.\n",BaseName()); 
    return 1;
  }
  if (dumpr_ && trajRadii_==NULL) {
    mprintf("Warning: setupTrajout[%s]: Radii are NULL and will not be printed.\n",
            BaseName());
    dumpr_=false;
  }
  return 0;
}

// PDBfile::SetWriteMode()
/** Set write mode to SINGLE, MODEL, or MULTI */
void PDBfile::SetWriteMode(PDBWRITEMODE modeIn) { 
  pdbWriteMode_ = modeIn; 
  //mprintf("PDB WRITE MODE SET TO %i\n",(int)pdbWriteMode_);
}

// PDBfile::SetDumpq()
void PDBfile::SetDumpq() {
  dumpq_ = true; 
  dumpr_ = true;
}

// PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int PDBfile::writeFrame(int set,double *X,double *V,double *box,double T) {
  char buffer[BUF_SIZE];
  int i,i3,res,atom,mol,lastAtomInMol,bufferSize;
  float Occ, B;

  // If writing 1 pdb per frame set up output filename and open
  if (pdbWriteMode_==MULTI) {
    NumberFilename(buffer,(char*)Name(),set + OUTPUTFRAMESHIFT);
    if (IO->Open(buffer,"wb")) return 1;
  // If specified, write MODEL keyword
  } else if (pdbWriteMode_==MODEL) {
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    Printf("MODEL     %i\n",set + OUTPUTFRAMESHIFT);
  }

  res=0; Occ=0.0; B=0.0;
  i3=0;
  atom=1; // Actual PDB atom number
  mol=0;
  lastAtomInMol=-1;
  if (trajAtomsPerMol_!=NULL)
    lastAtomInMol=trajAtomsPerMol_[0];
  for (i=0; i<pdbAtom_; i++) {
    // If this atom belongs to a new molecule print a TER card
    if (i == lastAtomInMol) {
      bufferSize=pdb_write_ATOM(buffer,PDBTER,atom,(char*)"",trajResNames_[res],
                                chainID_[i],res+1,0,0,0,0,0,(char*)"\0",dumpq_);
      IO->Write(buffer,bufferSize);
      atom++;
      mol++;
      lastAtomInMol += trajAtomsPerMol_[mol];
    }
    // figure out the residue number
    if ( i==trajResNums_[res+1] ) res++;
    if (dumpq_) Occ = (float) trajCharges_[i];
    if (dumpr_) B = (float) trajRadii_[i]; 
    bufferSize=pdb_write_ATOM(buffer,PDBATOM,atom,pdbAtomNames_[i],trajResNames_[res],
                              chainID_[i],res+1,X[i3],X[i3+1],X[i3+2],Occ,B,
                              (char*)"\0",dumpq_);
    IO->Write(buffer,bufferSize); 
    i3+=3;
    atom++;
  }

  // If writing 1 pdb per frame, close output file
  if (pdbWriteMode_==MULTI) {
    Printf("%-6s\n","END");
    IO->Close();
  // If MODEL keyword was written, write corresponding ENDMDL record
  } else if (pdbWriteMode_==MODEL) {
    Printf("ENDMDL\n");
  }

  return 0;
}

// PDBfile::info()
void PDBfile::info() {
  mprintf("is a PDB file");
  if (access_==WRITE) {
    if (pdbWriteMode_==MULTI)
      mprintf(" (1 file per frame)");
    else if (pdbWriteMode_==MODEL)
      mprintf(" (1 MODEL per frame)");
    if (dumpq_ && !dumpr_) 
      mprintf(", writing charges to occupancy column");
    else if (dumpr_ && !dumpq_) 
      mprintf(", writing GB radii to B-factor column");
    else if (dumpr_ && dumpq_)
      mprintf(", writing charges/GB radii to occupancy/B-factor columns");
  }
}
