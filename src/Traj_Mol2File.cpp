// Traj_Mol2File
#include "Traj_Mol2File.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_Mol2File::Traj_Mol2File() : 
  mol2WriteMode_(SINGLE),
  mol2Top_(0),
  hasCharges_(false)
{}

// Traj_Mol2File::openTraj()
int Traj_Mol2File::openTraj() {
  int err;

  err=0;
  switch (access_) {
    case READ : err = OpenFile(); break;
    case APPEND :
      mprinterr("Error: Append not supported for mol2 files.\n");
      err=1;
      break;
    case WRITE :
      // Set up title
      if (title_.empty())
        title_.assign("Cpptraj generated mol2 file.");
      // If writing 1 mol2 per frame do not open here
      if (mol2WriteMode_!=MULTI) err = OpenFile();
      break;
  }
      
  return err;
}

// Traj_Mol2File::closeTraj() {
void Traj_Mol2File::closeTraj() {
  // On WRITE only close if not writing 1 mol2 per frame
  if (access_==WRITE) {
    if (mol2WriteMode_!=MULTI) CloseFile();
  } else {
  // Otherwise just close the file
    CloseFile();
  }
}

// Traj_Mol2File::setupTrajin()
/** See how many MOLECULE records are in file, make sure num atoms match
  * parm and each frame.
  */
int Traj_Mol2File::setupTrajin(Topology *trajParm) {
  int frameAtom;
  int Frames=0;

  if (this->openTraj()) return -1;

  // Get @<TRIPOS>MOLECULE information for first frame
  if (ReadMolecule( IO )) return -1;
  // Check #atoms in mol2 file against #atoms in parm
  if (Mol2Natoms() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Mol2 file %s frame %i (%i) does not\n",
              BaseName(),Frames+OUTPUTFRAMESHIFT,Mol2Natoms());
    mprinterr("       match number in associated parmtop (%i)!\n",trajParm->Natom());
    return -1;
  }
  // Set title
  title_.assign( Mol2Title().c_str() );
  ++Frames;

  // See how many more MOLECULE records there are. Ensure same #atoms 
  // in each frame.
  while ( (frameAtom=NextMolecule( IO ))!=-1 ) {
    if (frameAtom != Mol2Natoms()) {
      mprintf("Warning: # atoms in Mol2 file %s frame %i (%i) not equal\n",
                BaseName(),Frames+OUTPUTFRAMESHIFT,frameAtom);
      mprintf("         to # atoms int first frame (%i).\n",Mol2Natoms());
      mprintf("         Only using frames 1-%i.\n",Frames);
      break;
    }
    ++Frames;
  }
 
  this->closeTraj();
  if (debug_>0) mprintf("      Mol2 file %s has %i frames.\n",BaseName(),Frames);
  seekable_ = false;
  hasTemperature_=false;
  hasBox_=false;

  return Frames;
}

// Traj_Mol2File::readFrame()
int Traj_Mol2File::readFrame(int set,double *X, double *V,double *box, double *T) {
  // Position file at next @<TRIPOS>ATOM tag
  if (ScanTo(IO, ATOM)) return 1;

  double *Xptr = X; 
  for (int atom = 0; atom < Mol2Natoms(); atom++) {
    if ( GetLine( IO ) ) return 1;
    Mol2XYZ(Xptr);
    //F->printAtomCoord(atom);
    Xptr += 3;
  }

  return 0;
}

// Traj_Mol2File::processWriteArgs()
int Traj_Mol2File::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("single")) this->SetWriteMode(MOL);
  if (argIn->hasKey("multi"))  this->SetWriteMode(MULTI);
  return 0;
}

// Traj_Mol2File::SetWriteMode()
/** Set write mode to SINGLE, MOL, or MULTI */
void Traj_Mol2File::SetWriteMode(MOL2WRITEMODE modeIn) {
  mol2WriteMode_ = modeIn;
  //mprintf("MOL2 WRITE MODE SET TO %i\n",(int)mol2WriteMode_);
}

// Traj_Mol2File::setupTrajout()
/** Set parm information required for write, and check write mode against
  * number of frames to be written.
  */
int Traj_Mol2File::setupTrajout(Topology *trajParm) {
  if (trajParm==NULL) return 1;
  mol2Top_ = trajParm;
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use @<TRIPOS>MOLECULE keyword to separate frames.
  if (mol2WriteMode_==SINGLE && mol2Top_->Nframes()>1) 
    mol2WriteMode_=MOL;

  // Set # atoms; if more than 99999 atoms the file may not write correctly
  SetMol2Natoms( mol2Top_->Natom() );
  if (Mol2Natoms() > 99999) {
    mprintf("Warning: %s: Large # of atoms (%i > 99999) for mol2 format.\n",
            BaseName(), Mol2Natoms());
    mprintf("         File may not write correctly.\n");
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

  // Setup output array for bonds. Bonds in parm are stored as:
  //   0_idx1 0_idx2 0_pidx 1_idx1 1_idx2 1_pidx ...
  // where idx is the atom index, i.e. atom# * 3, and pidx is the index
  // into the bond parm array. Create an array with atom #s only starting
  // from 1.
  trajBonds_.clear();
  for (Topology::bond_iterator bidx = mol2Top_->BondsStart();
                               bidx != mol2Top_->BondsEnd(); bidx++)
  {
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
  }
  for (Topology::bond_iterator bidx = mol2Top_->BondsH_Start(); 
                               bidx != mol2Top_->BondsH_End(); bidx++)
  {
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
    trajBonds_.push_back( ((*bidx)/3) + 1 );
    ++bidx;
  }
  // Check number of bonds
  if (trajBonds_.empty()) 
    mprintf("Warning: %s: topology does not contain bond information.\n",BaseName());
  else
    SetMol2Nbonds( trajBonds_.size() / 2 );

  return 0;
}

// Traj_Mol2File::writeFrame()
int Traj_Mol2File::writeFrame(int set, double *X, double *V,double *box, double T) {
  //mprintf("DEBUG: Calling Traj_Mol2File::writeFrame for set %i\n",set);
  if (mol2WriteMode_==MULTI) {
    std::string fname = NumberFilename(FullPathName(), set+OUTPUTFRAMESHIFT);
    if (IO->Open((char*)fname.c_str(),"wb")) return 1;
  }
  //@<TRIPOS>MOLECULE section
  // TODO: Add to Mol2File?
  Printf("@<TRIPOS>MOLECULE\n");
  // mol_name
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  // mol_type
  // charge_type
  // [status_bits
  // [mol_comment]]
  Printf("%s\n",title_.c_str());
  Printf("%5i %5i %5i %5i %5i\n",Mol2Natoms(), Mol2Nbonds(), 1, 0, 0);
  Printf("SMALL\n"); // May change this later
  if ( hasCharges_ )
    Printf("USER_CHARGES\n"); // May change this later
  else
    Printf("NO_CHARGES\n");
  Printf("\n\n");

  //@<TRIPOS>ATOM section
  Printf("@<TRIPOS>ATOM\n");
  double *Xptr = X;
  for (int i=0; i < Mol2Natoms(); i++) {
    const Atom atom = (*mol2Top_)[i];
    // figure out the residue number
    int res = atom.ResNum();
    // If atom type is blank, set to atom name.
    NameType atype = atom.Type();
    if ( atype == "" )
      atype = atom.Name();
    Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
                     i+1, atom.c_str(), Xptr[0], Xptr[1], Xptr[2],
                     *atype, res+1, mol2Top_->Res(res).c_str(), atom.Charge());
    Xptr += 3;
  }

  //@<TRIPOS>BOND section
  if (!trajBonds_.empty()) {
    Printf("@<TRIPOS>BOND\n");
    int bondnum = 1;
    for (std::vector<int>::iterator bond = trajBonds_.begin();
                                    bond != trajBonds_.end();
                                    bond++)
    {
      int firstBondAtom = *bond;
      ++bond;
      Printf("%5d %5d %5d 1\n",bondnum++,firstBondAtom, *bond);
    }
  }

  //@<TRIPOS>SUBSTRUCTURE section
  Printf("@<TRIPOS>SUBSTRUCTURE\n");
  int resnum = 1;
  for (Topology::res_iterator Res = mol2Top_->ResStart(); Res!=mol2Top_->ResEnd(); Res++) {
    Printf("%7d %4s %14d ****               0 ****  **** \n",
                      resnum++, (*Res).c_str(), (*Res).FirstAtom()+1);
  }

  // If writing 1 pdb per frame, close output file
  if (mol2WriteMode_==MULTI)
    IO->Close();

  return 0;
}
 
// Traj_Mol2File::info()
void Traj_Mol2File::info() {
  mprintf("is a Tripos Mol2 file");
  if (access_==WRITE) {
    if (mol2WriteMode_==MULTI)
      mprintf(" (1 file per frame)");
    else if (mol2WriteMode_==MOL)
      mprintf(" (1 MOLECULE per frame)");
  }
}
