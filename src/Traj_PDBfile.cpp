// Traj_PDBfile
#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_PDBfile::Traj_PDBfile() {
  pdbAtom_=0;
  pdbWriteMode_=SINGLE;
  dumpq_=false;
  dumpr_=false;
  chainchar_=' ';
}

//------------------------------------------------------------------------
// Traj_PDBfile::closeTraj()
void Traj_PDBfile::closeTraj() {
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

// Traj_PDBfile::openTraj()
int Traj_PDBfile::openTraj() {
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

// Traj_PDBfile::setupTrajin()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int Traj_PDBfile::setupTrajin(Topology *trajParm) {
  char buffer[BUF_SIZE_];
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
  if (debug_>0) mprintf("Traj_PDBfile: %s has %i atoms, %i frames.\n",BaseName(),
                       pdbAtom_,Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",BaseName(),
            numMismatch,trajParm->parmName);

  return Frames;
}

// Traj_PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int Traj_PDBfile::readFrame(int set,double *X, double *V,double *box, double *T) {
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

// Traj_PDBfile::processWriteArgs()
int Traj_PDBfile::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("dumpq")) this->SetDumpq();
  if (argIn->hasKey("model")) this->SetWriteMode(MODEL);
  if (argIn->hasKey("multi")) this->SetWriteMode(MULTI);
  char *temp = argIn->getKeyString("chainid",NULL);
  if (temp!=NULL) chainchar_ = temp[0];
  return 0;
}

// Traj_PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int Traj_PDBfile::setupTrajout(Topology *trajParm) {
  pdbTop = trajParm;

  pdbAtom_ = pdbTop->Natom();

  // Set a chainID for each atom
  chainID_.clear();
  chainID_.resize(pdbAtom_, chainchar_);

  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (pdbWriteMode_==SINGLE && trajParm->Nframes()>1) 
    pdbWriteMode_=MODEL;
  return 0;
}

// Traj_PDBfile::SetWriteMode()
/** Set write mode to SINGLE, MODEL, or MULTI */
void Traj_PDBfile::SetWriteMode(PDBWRITEMODE modeIn) { 
  pdbWriteMode_ = modeIn; 
  //mprintf("PDB WRITE MODE SET TO %i\n",(int)pdbWriteMode_);
}

// Traj_PDBfile::SetDumpq()
void Traj_PDBfile::SetDumpq() {
  dumpq_ = true; 
  dumpr_ = true;
}

// Traj_PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int Traj_PDBfile::writeFrame(int set,double *X,double *V,double *box,double T) {
  char fname[BUF_SIZE_]; // TODO: rewrite with string or sstream
  float Occ, B;

  // If writing 1 pdb per frame set up output filename and open
  if (pdbWriteMode_==MULTI) {
    NumberFilename(fname,(char*)Name(),set + OUTPUTFRAMESHIFT);
    if (IO->Open(buffer,"wb")) return 1;
  // If specified, write MODEL keyword
  } else if (pdbWriteMode_==MODEL) {
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    Printf("MODEL     %i\n",set + OUTPUTFRAMESHIFT);
  }

  Occ = 0.0; 
  B = 0.0;
  int anum = 1; // Actual PDB atom number
  Topology::mol_iterator mol = pdbParm->MolStart();
  int lastAtomInMol = (*mol).EndAtom();
  double *Xptr = X;
  for (int i = 0; i < pdbAtom_; i++) {
    const Atom atom = (*pdbParm_)[i];
    int res = atom.ResNum();
    // If this atom belongs to a new molecule print a TER card
    if (i == lastAtomInMol) {
      pdb_write_ATOM(IO, PDBTER, anum, "", pdbTop->Res(res).Name(),
                     chainID_[i], res+1, 0, 0, 0, 0, 0, (char*)"\0", dumpq_);
      ++anum;
      ++mol;
      lastAtomInMol = (*mol).EndAtom();
    }
    if (dumpq_) Occ = (float) atom.Charge();
    if (dumpr_) B = (float) atom.Radius();
    pdb_write_ATOM(IO, PDBATOM, atom, atom.Name(), pdbTop->Res(res).Name(),
                   chainID_[i], res+1, Xptr[0], Xptr[1], Xptr[2], Occ, B, (char*)"\0", dumpq_);
    Xptr += 3;
    ++atom;
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

// Traj_PDBfile::info()
void Traj_PDBfile::info() {
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
