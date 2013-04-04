// Traj_PDBfile
#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// CONSTRUCTOR
Traj_PDBfile::Traj_PDBfile() :
  pdbAtom_(0),
  ter_num_(0),
  pdbWriteMode_(NONE),
  dumpq_(false),
  dumpr_(false),
  pdbTop_(0),
  chainchar_(' ')
{}

//------------------------------------------------------------------------
bool Traj_PDBfile::ID_TrajFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}

// Traj_PDBfile::closeTraj()
/** If not writing one PDB per frame write the END record. */
void Traj_PDBfile::closeTraj() {
  if ( (pdbWriteMode_ == SINGLE || pdbWriteMode_ == MODEL) &&
        file_.IsOpen() )
    file_.Printf("%-6s\n","END");
  if (pdbWriteMode_ != MULTI)
    file_.CloseFile();
}

// Traj_PDBfile::openTrajin()
int Traj_PDBfile::openTrajin() {
  return file_.OpenFile();
}

// Traj_PDBfile::setupTrajin()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int Traj_PDBfile::setupTrajin(std::string const& fname, Topology* trajParm)
{
  int atom;
  pdbWriteMode_ = NONE;
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  // Two strats - check for MODEL keywords or see how many times natom ATOM
  // records can be read. Currently employing the latter.
  int Frames = 0;
  int numMismatch = 0;
  bool scanPDB = true;
  while (scanPDB) {
    atom = 0;
    while ( atom < trajParm->Natom() ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( file_.NextLine() == 0 ) {
        scanPDB = false;
        break;
      }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (!file_.IsPDBatomKeyword()) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        Atom pdbatom = file_.pdb_Atom();
        if ( pdbatom.Name() != (*trajParm)[atom].Name() ) {
          if (debug_>1) 
            mprintf("Warning: %s: PDB atom %i name [%s] does not match parm atom name [%s]\n",
                    file_.Filename().base(), atom+1, *(pdbatom.Name()), 
                    *((*trajParm)[atom].Name()));
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
        mprintf("Warning: PDB %s: Reading frame %i, got %i atoms, expected %i.\n",
                file_.Filename().base(), Frames+1, atom, pdbAtom_);
        mprintf("Warning: Only using frames 1-%i\n", Frames);
        scanPDB = false;
        break;
      }
    }  
    if (scanPDB) ++Frames;
  }
  file_.CloseFile(); 
  if (Frames<1) {
    mprinterr("Error: PDB %s: No frames read. atom=%i expected %i.\n",
              file_.Filename().base(), atom, trajParm->Natom());
    return TRAJIN_ERR;
  }
  if (debug_>0) mprintf("Traj_PDBfile: %s has %i atoms, %i frames.\n",
                        file_.Filename().base(), pdbAtom_, Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",
            file_.Filename().base(), numMismatch, trajParm->c_str());
  return Frames;
}

// Traj_PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int Traj_PDBfile::readFrame(int set,double *X, double *V,double *box, double *T) 
{
  int atom = 0;
  double *Xptr = X; 
  while (atom < pdbAtom_) {
    if ( file_.NextLine() == 0 ) return 1;
    // Skip non-ATOM records
    if (file_.IsPDBatomKeyword()) {
      // Read current PDB record XYZ into Frame
      file_.pdb_XYZ( Xptr );
      ++atom; 
      Xptr += 3;
    }
  }

  return 0;
}

// Traj_PDBfile::processWriteArgs()
int Traj_PDBfile::processWriteArgs(ArgList& argIn) {
  pdbWriteMode_ = SINGLE;
  if (argIn.hasKey("dumpq")) {
   dumpq_ = true; 
   dumpr_ = true;
  }
  if (argIn.hasKey("teradvance")) ter_num_ = 1;
  if (argIn.hasKey("model")) pdbWriteMode_ = MODEL;
  if (argIn.hasKey("multi")) pdbWriteMode_ = MULTI;
  std::string temp = argIn.GetStringKey("chainid");
  if (!temp.empty()) chainchar_ = temp[0];
  return 0;
}

// Traj_PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int Traj_PDBfile::setupTrajout(std::string const& fname, Topology* trajParm,
                               int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  pdbTop_ = trajParm;
  pdbAtom_ = pdbTop_->Natom();
  // Set up file
  if (append && pdbWriteMode_ != MULTI) {
    if (file_.SetupAppend( fname, debug_)) return 1;
  } else {
    if (append && pdbWriteMode_ == MULTI)
      mprintf("Warning: 'append' not compatible with 'multi' pdb write.\n");
    if (file_.SetupWrite( fname, debug_ )) return 1;
  }
  // Set a chainID for each atom
  // TODO: Set different chain ID for solute mols and solvent
  chainID_.clear();
  chainID_.resize(pdbAtom_, chainchar_);
  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (append || (pdbWriteMode_==SINGLE && NframesToWrite>1)) 
    pdbWriteMode_ = MODEL;
  // TODO: Setup title
  // Open here if writing to single file
  if (pdbWriteMode_ != MULTI)
    return file_.OpenFile();
  return 0;
}

// Traj_PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int Traj_PDBfile::writeFrame(int set,double *X,double *V,double *box,double T) {
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame set up output filename and open
    if (file_.OpenWriteWithName( NumberFilename(file_.Filename().Full(), set+1) )) return 1;
  } else if (pdbWriteMode_==MODEL) {
    // If specified, write MODEL keyword
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    file_.Printf("MODEL     %i\n", set + 1);
  }

  float Occ = 1.0; 
  float B = 0.0;
  int anum = 1; // Actual PDB atom number
  int aidx = 0; // Atom index in topology
  Topology::mol_iterator mol = pdbTop_->MolStart();
  int lastAtomInMol = (*mol).EndAtom();
  double *Xptr = X;
  for (Topology::atom_iterator atom = pdbTop_->begin(); atom != pdbTop_->end(); ++atom) {
    int res = (*atom).ResNum();
    // If this atom belongs to a new molecule print a TER card
    // Use res instead of res+1 since this TER belongs to last mol/res
    if (aidx == lastAtomInMol) {
      file_.WriteTER( anum, pdbTop_->Res(res-1).Name(), chainID_[aidx], res );
      anum += ter_num_;
      ++mol;
      lastAtomInMol = (*mol).EndAtom();
    }
    if (dumpq_) Occ = (float) (*atom).Charge();
    if (dumpr_) B = (float) (*atom).Radius();
    file_.WriteRec(PDBfile::ATOM, anum++, (*atom).Name(), pdbTop_->Res(res).Name(),
                   chainID_[aidx++], res+1, Xptr[0], Xptr[1], Xptr[2], Occ, B, 
                   (*atom).ElementName(), 0, dumpq_);
    Xptr += 3;
  }
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame, close output file
    file_.Printf("%-6s\n","END");
    file_.CloseFile();
  } else if (pdbWriteMode_==MODEL) {
    // If MODEL keyword was written, write corresponding ENDMDL record
    file_.Printf("ENDMDL\n");
  }

  return 0;
}

// Traj_PDBfile::Info()
void Traj_PDBfile::Info() {
  mprintf("is a PDB file");
  if (pdbWriteMode_ != NONE) {
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
