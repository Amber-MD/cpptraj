// Traj_Mol2File
#include "Traj_Mol2File.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// CONSTRUCTOR
Traj_Mol2File::Traj_Mol2File() : 
  mol2WriteMode_(NONE),
  mol2Top_(0),
  hasCharges_(false)
{}

bool Traj_Mol2File::ID_TrajFormat(CpptrajFile& fileIn) {
  return Mol2File::ID_Mol2( fileIn );
}

// Traj_Mol2File::openTrajin()
int Traj_Mol2File::openTrajin() {
  return file_.OpenFile(); 
}

// Traj_Mol2File::closeTraj() {
void Traj_Mol2File::closeTraj() {
  if ( mol2WriteMode_ != MULTI )
    file_.CloseFile();
}

// Traj_Mol2File::setupTrajin()
/** See how many MOLECULE records are in file, make sure num atoms match
  * parm and each frame.
  */
int Traj_Mol2File::setupTrajin(std::string const& fname, Topology* trajParm)
{
  int Frames=0;
  mol2WriteMode_ = NONE;
  if (file_.SetupRead( fname, debug_)) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  // Get @<TRIPOS>MOLECULE information for first frame
  if (file_.ReadMolecule()) return TRAJIN_ERR;
  // Check #atoms in mol2 file against #atoms in parm
  if (file_.Mol2Natoms() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Mol2 file %s frame %i (%i) does not\n",
              file_.Filename().base(), Frames+1, file_.Mol2Natoms());
    mprinterr("Error: match number in associated parmtop (%i)!\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Set title
  SetTitle( file_.Mol2Title() );
  Frames = 1;
  // See how many more MOLECULE records there are. Ensure same #atoms 
  // in each frame.
  int frameAtom;
  while ( (frameAtom = file_.NextMolecule( ))!=-1 ) {
    if (frameAtom != file_.Mol2Natoms()) {
      mprintf("Warning: # atoms in Mol2 file %s frame %i (%i) not equal\n",
              file_.Filename().base(), Frames+1, frameAtom);
      mprintf("Warning: to # atoms int first frame (%i).\n", file_.Mol2Natoms());
      mprintf("Warning: Only using frames 1-%i.\n", Frames);
      break;
    }
    ++Frames;
  }
  file_.CloseFile();
  if (debug_>0) mprintf("\tMol2 file %s has %i frames.\n",file_.Filename().base(), Frames);

  return Frames;
}

// Traj_Mol2File::readFrame()
int Traj_Mol2File::readFrame(int set,double *X, double *V,double *box, double *T) {
  // Position file at next @<TRIPOS>ATOM tag
  if (file_.ScanTo(Mol2File::ATOM)) return 1;
  double *Xptr = X; 
  for (int atom = 0; atom < file_.Mol2Natoms(); atom++) {
    if (file_.Mol2XYZ(Xptr)) return 1;
    //F->printAtomCoord(atom);
    Xptr += 3;
  }
  return 0;
}

// Traj_Mol2File::processWriteArgs()
int Traj_Mol2File::processWriteArgs(ArgList& argIn) {
  mol2WriteMode_ = SINGLE; // Default to single writes
  if (argIn.hasKey("single")) mol2WriteMode_ = MOL;
  if (argIn.hasKey("multi"))  mol2WriteMode_ = MULTI;
  return 0;
}

// Traj_Mol2File::setupTrajout()
/** Set parm information required for write, and check write mode against
  * number of frames to be written.
  */
int Traj_Mol2File::setupTrajout(std::string const& fname, Topology* trajParm,
                                int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  mol2Top_ = trajParm;
  // Set up file
  if (append && mol2WriteMode_ != MULTI) {
    if (file_.SetupAppend( fname, debug_)) return 1;
  } else {
    if (mol2WriteMode_ == MULTI)
      mprintf("Warning: 'append' not compatible with 'multi' mol2 write.\n");
    if (file_.SetupWrite( fname, debug_ )) return 1;
  }
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use @<TRIPOS>MOLECULE keyword to separate frames.
  if (append || (mol2WriteMode_==SINGLE && NframesToWrite>1)) 
    mol2WriteMode_ = MOL;
  // Set # atoms; if more than 99999 atoms the file may not write correctly
  file_.SetMol2Natoms( mol2Top_->Natom() );
  if (file_.Mol2Natoms() > 99999) {
    mprintf("Warning: %s: Large # of atoms (%i > 99999) for mol2 format.\n",
            file_.Filename().base(), file_.Mol2Natoms());
    mprintf("Warning: File may not write correctly.\n");
  }
  // TODO: Change this, right now for backwards compat. only!
  // If all charges == 0 set noCharges.
  hasCharges_ = false;
  for (Topology::atom_iterator atom = mol2Top_->begin(); atom != mol2Top_->end(); atom++)
  {
    if ( (*atom).Charge() != 0 ) {
      hasCharges_ = true;
      break;
    }
  }
  // Set Title
  if (Title().empty())
    SetTitle("Cpptraj generated mol2 file.");
  file_.SetMol2Title( Title() );
  // Setup output array for bonds. Bonds in parm are stored as:
  //   0_idx1 0_idx2 0_pidx 1_idx1 1_idx2 1_pidx ...
  // where idx is the atom index, i.e. atom# * 3, and pidx is the index
  // into the bond parm array. Create an array with atom #s only starting
  // from 1.
  trajBonds_.clear();
  for (std::vector<int>::const_iterator bidx = mol2Top_->Bonds().begin();
                                        bidx != mol2Top_->Bonds().end(); bidx++)
  {
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
  }
  for (std::vector<int>::const_iterator bidx = mol2Top_->BondsH().begin();
                                        bidx != mol2Top_->BondsH().end(); bidx++)
  {
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
  }
  // Check number of bonds
  if (trajBonds_.empty()) 
    mprintf("Warning: %s: topology does not contain bond information.\n",
            file_.Filename().base());
  else
    file_.SetMol2Nbonds( trajBonds_.size() / 2 );
  // Open here if writing to a single file
  if ( mol2WriteMode_ != MULTI )
    return file_.OpenFile(); 
  return 0;
}

// Traj_Mol2File::writeFrame()
int Traj_Mol2File::writeFrame(int set, double *X, double *V,double *box, double T) {
  //mprintf("DEBUG: Calling Traj_Mol2File::writeFrame for set %i\n",set);
  if (mol2WriteMode_==MULTI) {
    if (file_.OpenWriteWithName( NumberFilename( file_.Filename().Full(), set+1) )) return 1;
  }
  //@<TRIPOS>MOLECULE section
  file_.WriteMolecule( hasCharges_, mol2Top_->Nres() );
  //@<TRIPOS>ATOM section
  file_.WriteHeader(Mol2File::ATOM);
  double *Xptr = X;
  int atnum = 1;
  for (Topology::atom_iterator atom = mol2Top_->begin(); atom != mol2Top_->end(); ++atom) {
    // figure out the residue number
    int res = (*atom).ResNum();
    // If atom type is blank, set to atom name.
    NameType atype = (*atom).Type();
    if ( atype == "" )
      atype = (*atom).Name();
    file_.Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
                     atnum++, (*atom).c_str(), Xptr[0], Xptr[1], Xptr[2],
                     *atype, res+1, mol2Top_->Res(res).c_str(), (*atom).Charge());
    Xptr += 3;
  }
  //@<TRIPOS>BOND section
  if (!trajBonds_.empty()) {
    file_.WriteHeader(Mol2File::BOND); 
    int bondnum = 1;
    for (std::vector<int>::iterator bond = trajBonds_.begin();
                                    bond != trajBonds_.end();
                                    bond++)
    {
      int firstBondAtom = *bond;
      ++bond;
      file_.Printf("%5d %5d %5d 1\n",bondnum++,firstBondAtom, *bond);
    }
  }
  //@<TRIPOS>SUBSTRUCTURE section
  file_.WriteHeader(Mol2File::SUBSTRUCT);
  int resnum = 1;
  for (Topology::res_iterator Res = mol2Top_->ResStart(); Res!=mol2Top_->ResEnd(); Res++) {
    file_.Printf("%7d %4s %14d ****               0 ****  **** \n",
                      resnum++, (*Res).c_str(), (*Res).FirstAtom()+1);
  }
  // If writing 1 mol2 per frame, close output file
  if (mol2WriteMode_==MULTI)
    file_.CloseFile();

  return 0;
}
 
// Traj_Mol2File::info()
void Traj_Mol2File::Info() {
  mprintf("is a Tripos Mol2 file");
  if (mol2WriteMode_==MULTI)
    mprintf(" (1 file per frame)");
  else if (mol2WriteMode_==MOL)
    mprintf(" (1 MOLECULE per frame)");
}
