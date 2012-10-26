// Traj_PDBfile
#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// CONSTRUCTOR
Traj_PDBfile::Traj_PDBfile() :
  pdbAtom_(0),
  pdbWriteMode_(SINGLE),
  dumpq_(false),
  dumpr_(false),
  pdbTop_(0),
  chainchar_(' ')
{}

//------------------------------------------------------------------------
bool Traj_PDBfile::ID_TrajFormat() {
  // Assumes already set up
  if (OpenFile()) return false;
  bool ispdbfile = ID( IO );
  CloseFile();
  return ispdbfile;
}

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
      mprinterr("Error: %s: Append currently not supported for PDB files.\n",BaseFileStr());
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
  int atom;
  if ( this->openTraj() ) return -1;

  // Two strats - check for MODEL keywords or see how many times natom ATOM
  // records can be read.
  // Currently employing the latter.
  int Frames = 0;
  int numMismatch = 0;
  bool scanPDB = true;
  while (scanPDB) {
    atom = 0;
    while ( atom < trajParm->Natom() ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( !PDB_GetNextRecord( IO ) ) {
        scanPDB = false;
        break;
      }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (!IsPDBatomKeyword()) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        Atom pdbatom = pdb_Atom();
        if ( pdbatom.Name() != (*trajParm)[atom].Name() ) {
          if (debug_>1) 
            mprintf("Warning: %s: PDB atom %i name [%s] does not match parm atom name [%s]\n",
                    BaseFileStr(),atom+1,*(pdbatom.Name()),*((*trajParm)[atom].Name()));
          ++numMismatch;
        }
      }
      ++atom;
    }
    if (Frames==0) {
      // First frame #atoms 
      pdbAtom_ = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom_!=atom) {
        mprintf("Warning: PDB %s: Reading frame %i, got %i atoms, expected %i.\n",BaseFileStr(),
                  Frames+OUTPUTFRAMESHIFT,atom,pdbAtom_);
        mprintf("         Only using frames 1-%i\n",Frames);
        scanPDB = false;
        break;
      }
    }  
    if (scanPDB) ++Frames;
  }
  this->closeTraj();

  if (Frames<1) {
    mprinterr("Error: PDB %s: No frames read. atom=%i expected %i.\n",BaseFileStr(),
            atom,trajParm->Natom());
    return -1;
  }
  if (debug_>0) mprintf("Traj_PDBfile: %s has %i atoms, %i frames.\n",BaseFileStr(),
                       pdbAtom_,Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",BaseFileStr(),
            numMismatch,trajParm->c_str());

  return Frames;
}

// Traj_PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int Traj_PDBfile::readFrame(int set,double *X, double *V,double *box, double *T) {
  int atom = 0;
  double *Xptr = X; 
  while (atom < pdbAtom_) {
    if ( !PDB_GetNextRecord( IO ) ) return 1;
    // Skip non-ATOM records
    if (!IsPDBatomKeyword()) continue;
    // Read current PDB record XYZ into Frame
    pdb_XYZ( Xptr );
    ++atom; 
    Xptr += 3;
  }

  return 0;
}

// Traj_PDBfile::processWriteArgs()
int Traj_PDBfile::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("dumpq")) this->SetDumpq();
  if (argIn->hasKey("model")) this->SetWriteMode(MODEL);
  if (argIn->hasKey("multi")) this->SetWriteMode(MULTI);
  ArgList::ConstArg temp = argIn->getKeyString("chainid");
  if (temp!=NULL) chainchar_ = temp[0];
  return 0;
}

// Traj_PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int Traj_PDBfile::setupTrajout(Topology *trajParm, int NframesToWrite) {
  pdbTop_ = trajParm;

  pdbAtom_ = pdbTop_->Natom();

  // Set a chainID for each atom
  chainID_.clear();
  chainID_.resize(pdbAtom_, chainchar_);

  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (pdbWriteMode_==SINGLE && NframesToWrite>1) 
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
  // If writing 1 pdb per frame set up output filename and open
  if (pdbWriteMode_==MULTI) {
    std::string fname = NumberFilename(FullFileName(), set + OUTPUTFRAMESHIFT);
    //mprintf("DEBUG: Traj_PDBfile::writeFrame: set %i, MULTI [%s]\n",set,fname.c_str());
    if (IO->Open(fname.c_str(), "wb")) return 1;
  // If specified, write MODEL keyword
  } else if (pdbWriteMode_==MODEL) {
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    Printf("MODEL     %i\n",set + OUTPUTFRAMESHIFT);
  }

  float Occ = 0.0; 
  float B = 0.0;
  int anum = 1; // Actual PDB atom number
  Topology::mol_iterator mol = pdbTop_->MolStart();
  int lastAtomInMol = (*mol).EndAtom();
  double *Xptr = X;
  for (int i = 0; i < pdbAtom_; i++) {
    const Atom atom = (*pdbTop_)[i];
    int res = atom.ResNum();
    // If this atom belongs to a new molecule print a TER card
    // Use res instead of res+1 since this TER belongs to last mol/res
    if (i == lastAtomInMol) {
      pdb_write_ATOM(IO, PDBTER, anum, "", pdbTop_->Res(res-1).Name(),
                     chainID_[i], res, 0, 0, 0, 0, 0, "", dumpq_);
      ++anum;
      ++mol;
      lastAtomInMol = (*mol).EndAtom();
    }
    if (dumpq_) Occ = (float) atom.Charge();
    if (dumpr_) B = (float) atom.Radius();
    pdb_write_ATOM(IO, PDBATOM, anum, atom.Name(), pdbTop_->Res(res).Name(),
                   chainID_[i], res+1, Xptr[0], Xptr[1], Xptr[2], Occ, B, "", dumpq_);
    Xptr += 3;
    ++anum;
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
