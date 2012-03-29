// Mol2File
#include <cstdio> // sscanf
#include "Traj_Mol2File.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Mol2File::Mol2File() {
  mol2atom_=0;
  mol2bonds_=0;
  mol2WriteMode_=SINGLE;

  trajnres_ = 0;
  trajAtomNames_ = NULL; 
  trajTypes_ = NULL;
  trajResNames_ = NULL; 
  trajResNums_ = NULL;
  trajCharges_ = NULL;
}

// Mol2File::openTraj()
int Mol2File::openTraj() {
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

// Mol2File::closeTraj() {
void Mol2File::closeTraj() {
  // On WRITE only close if not writing 1 mol2 per frame
  if (access_==WRITE) {
    if (mol2WriteMode_!=MULTI) CloseFile();
  } else {
  // Otherwise just close the file
    CloseFile();
  }
}

// Mol2File::setupTrajin()
/** See how many MOLECULE records are in file, make sure num atoms match
  * parm and each frame.
  */
int Mol2File::setupTrajin(AmberParm *trajParm) {
  int frameAtom;
  int Frames=0;
  char buffer[MOL2BUFFERSIZE];

  if (this->openTraj()) return -1;

  // Get @<TRIPOS>MOLECULE information
  while (Mol2ScanTo(IO, MOLECULE)==0) {
    //   Scan title
    if ( IO->Gets(buffer,MOL2BUFFERSIZE) ) return -1; 
    // On first frame set title
    if (Frames==0) {
      title_.assign(buffer);
      if (debug_>0) mprintf("    Mol2 Title: [%s]\n",title_.c_str());
    }
    //   Scan # atoms
    if ( IO->Gets(buffer,MOL2BUFFERSIZE) ) return -1;
    sscanf(buffer,"%i",&frameAtom);
    if (Frames==0) {
      mol2atom_=frameAtom;
      if (debug_>0) mprintf("    Mol2 #atoms: %i\n",mol2atom_);
    }
    if (Frames==0 && frameAtom!=trajParm->natom) {
      mprinterr("Error: Number of atoms in Mol2 file %s frame %i (%i) does not\n",
              BaseName(),Frames+OUTPUTFRAMESHIFT,frameAtom);
      mprinterr("       match number in associated parmtop (%i)!\n",trajParm->natom);
      return -1;
    }
    if (frameAtom!=mol2atom_) {
      mprintf("Warning: # atoms in Mol2 file %s frame %i (%i) not equal\n",
                BaseName(),Frames+OUTPUTFRAMESHIFT,frameAtom);
      mprintf("         to # atoms int first frame (%i).\n",mol2atom_);
      mprintf("         Only using frames 1-%i.\n",Frames);
      break;
    }
    Frames++;
  }

  this->closeTraj();
  if (debug_>0) mprintf("      Mol2 file %s has %i frames.\n",BaseName(),Frames);
  seekable_ = false;
  hasTemperature_=false;
  hasBox_=false;

  return Frames;
}

// Mol2File::readFrame()
int Mol2File::readFrame(int set,double *X, double *V,double *box, double *T) {
  char buffer[MOL2BUFFERSIZE];
  int atom, atom3;

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(IO, ATOM)) return 1;

  atom3=0;
  for (atom = 0; atom < mol2atom_; atom++) {
    if (IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    sscanf(buffer,"%*i %*s %lf %lf %lf",X+atom3, X+atom3+1, X+atom3+2);
    //F->printAtomCoord(atom);
    atom3+=3;
  }

  return 0;
}

// Mol2File::processWriteArgs()
int Mol2File::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("single")) this->SetWriteMode(MOL);
  if (argIn->hasKey("multi"))  this->SetWriteMode(MULTI);
  return 0;
}

// Mol2File::SetWriteMode()
/** Set write mode to SINGLE, MOL, or MULTI */
void Mol2File::SetWriteMode(MOL2WRITEMODE modeIn) {
  mol2WriteMode_ = modeIn;
  //mprintf("MOL2 WRITE MODE SET TO %i\n",(int)mol2WriteMode_);
}

// Mol2File::setupTrajout()
/** Set parm information required for write, and check write mode against
  * number of frames to be written.
  */
int Mol2File::setupTrajout(AmberParm *trajParm) {
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use @<TRIPOS>MOLECULE keyword to separate frames.
  if (mol2WriteMode_==SINGLE && trajParm->parmFrames>1) 
    mol2WriteMode_=MOL;

  // Set # atoms; if more than 99999 atoms the file may not write correctly
  mol2atom_ = trajParm->natom;
  if (mol2atom_>99999) {
    mprintf("Warning: %s: Large # of atoms (%i > 99999) for mol2 format.\n",
            BaseName(),mol2atom_);
    mprintf("         File may not write correctly.\n");
  }

  // Check number of bonds
  if ( trajParm->BondArray( trajBonds_ ) ) 
    mprintf("Warning: %s: topology does not contain bond information.\n",BaseName());
  // The BondArray function results in an array consisting of all bonds
  // with both atoms of the bond, hence the number of bonds is size / 2
  mol2bonds_ = (int)trajBonds_.size();
  mol2bonds_ /= 2;

  // Set information from parm
  trajnres_ = trajParm->Nres();
  trajAtomNames_ = trajParm->AtomNames_ptr();
  trajTypes_ = trajParm->AtomTypes_ptr();
  // If types are not set just use atom names
  if (trajTypes_ == NULL) trajTypes_ = trajAtomNames_;
  trajResNames_ = trajParm->ResidueNames_ptr();
  trajResNums_ = trajParm->ResAtomNums_ptr();
  trajCharges_ = trajParm->Charges_ptr();

  // Check that all parm info is indeed present
  if (trajAtomNames_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Atom names are NULL.\n",BaseName());
    return 1;
  }
  if (trajTypes_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Atom types are NULL.\n",BaseName());
    return 1;
  }
  if (trajResNames_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Residue names are NULL.\n",BaseName());
    return 1;
  }
  if (trajResNums_==NULL) {
    mprinterr("Error: setupTrajout [%s]: Residue #s are NULL.\n",BaseName());
    return 1;
  }
  //if (trajCharges_==NULL) {
  //  mprinterr("Error: setupTrajout [%s]: Charges are NULL.\n",BaseName());
  //  return 1;
  //}
  if (trajBonds_.empty()) {
    mprintf("Warning: setupTrajout [%s]: No bond information present in parm %s\n",
            BaseName(),trajParm->parmName);
  }
  return 0;
}

// Mol2File::writeFrame()
int Mol2File::writeFrame(int set, double *X, double *V,double *box, double T) {
  char buffer[1024];
  int atom, atom3, res;
  double Charge;

  //mprintf("DEBUG: Calling Mol2File::writeFrame for set %i\n",set);
  if (mol2WriteMode_==MULTI) {
    NumberFilename(buffer,(char*)Name(),set+OUTPUTFRAMESHIFT);
    if (IO->Open(buffer,"wb")) return 1;
  }
  //@<TRIPOS>MOLECULE section
  Printf("@<TRIPOS>MOLECULE\n");
  // mol_name
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  // mol_type
  // charge_type
  // [status_bits
  // [mol_comment]]
  Printf("%s\n",title_.c_str());
  Printf("%5i %5i %5i %5i %5i\n",mol2atom_,mol2bonds_,1,0,0);
  Printf("SMALL\n"); // May change this later
  if (trajCharges_!=NULL)
    Printf("USER_CHARGES\n"); // May change this later
  else
    Printf("NO_CHARGES\n");
  Printf("\n\n");

  //@<TRIPOS>ATOM section
  Charge = 0.0;
  Printf("@<TRIPOS>ATOM\n");
  atom3=0;
  res = 0;
  for (atom=0; atom < mol2atom_; atom++) {
    // figure out the residue number
    if ( atom==trajResNums_[res+1] ) res++;
    if (trajCharges_!=NULL) Charge = trajCharges_[atom]; 
    Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
                     atom+1, trajAtomNames_[atom], X[atom3], X[atom3+1], X[atom3+2],
                     trajTypes_[atom], res+1, trajResNames_[res], Charge);
    atom3+=3;
  }

  //@<TRIPOS>BOND section
  if (!trajBonds_.empty()) {
    // Atom #s in the bonds and bondh array are * 3
    Printf("@<TRIPOS>BOND\n");
    atom=1;
    for (std::vector<int>::iterator bond = trajBonds_.begin();
                                    bond != trajBonds_.end();
                                    bond++)
    {
      int firstBondAtom = *bond;
      ++bond;
      Printf("%5d %5d %5d 1\n",atom++,firstBondAtom + 1,*bond + 1);
    }
  }

  //@<TRIPOS>SUBSTRUCTURE section
  Printf("@<TRIPOS>SUBSTRUCTURE\n");
  for (res = 0; res < trajnres_; res++) {
    Printf("%7d %4s %14d ****               0 ****  **** \n",
                      res+1, trajResNames_[res],trajResNums_[res]+1);
  }

  // If writing 1 pdb per frame, close output file
  if (mol2WriteMode_==MULTI)
    IO->Close();

  return 0;
}
 
// Mol2File::info()
void Mol2File::info() {
  mprintf("is a Tripos Mol2 file");
  if (access_==WRITE) {
    if (mol2WriteMode_==MULTI)
      mprintf(" (1 file per frame)");
    else if (mol2WriteMode_==MOL)
      mprintf(" (1 MOLECULE per frame)");
  }
}
