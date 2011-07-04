// Mol2File
#include <cstdio> // sscanf
#include "Traj_Mol2File.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Mol2File::Mol2File() {
  mol2atom=0;
  mol2bonds=0;
  mol2WriteMode=SINGLE;

  trajnres = 0;
  trajnbondsh = 0;
  trajnbonds = 0;
  trajAtomNames = NULL; 
  trajTypes = NULL;
  trajResNames = NULL; 
  trajResNums = NULL;
  trajCharges = NULL;
  trajBonds = NULL;
  trajBondsh = NULL;
}

// DESTRUCTOR
Mol2File::~Mol2File() {
}

/* Mol2File::openTraj()
 */
int Mol2File::openTraj() {
  int err;

  err=0;
  switch (tfile->access) {
    case READ : err = tfile->OpenFile(); break;
    case APPEND :
      mprintf("Error: Append not supported for mol2 files.\n");
      err=1;
      break;
    case WRITE :
      // Set up title
      if (title==NULL)
        this->SetTitle((char*)"Cpptraj generated mol2 file.\0");
      // If writing 1 mol2 per frame do not open here
      if (mol2WriteMode!=MULTI) err = tfile->OpenFile();
      break;
  }
      
  return err;
}

/* Mol2File::closeTraj() {
 */
void Mol2File::closeTraj() {
  // On WRITE only close if not writing 1 mol2 per frame
  if (tfile->access==WRITE) {
    if (mol2WriteMode!=MULTI) tfile->CloseFile();
  } else {
  // Otherwise just close the file
    tfile->CloseFile();
  }
}

/* Mol2File::setupRead()
 * See how many MOLECULE records are in file, make sure num atoms match
 * parm and each frame.
 */
int Mol2File::setupRead(int natom) {
  int frameAtom;
  int Frames=0;
  char buffer[MOL2BUFFERSIZE];

  if (this->openTraj()) return -1;

  // Get @<TRIPOS>MOLECULE information
  while (Mol2ScanTo(this->tfile, MOLECULE)==0) {
    //   Scan title
    if ( tfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return -1; 
    // On first frame set title
    if (Frames==0) {
      this->SetTitle(buffer);
      if (debug>0) mprintf("    Mol2 Title: [%s]\n",title);
    }
    //   Scan # atoms
    if ( tfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return -1;
    sscanf(buffer,"%i",&frameAtom);
    if (Frames==0) {
      mol2atom=frameAtom;
      if (debug>0) mprintf("    Mol2 #atoms: %i\n",mol2atom);
    }
    if (frameAtom!=natom) {
      mprinterr("Error: Number of atoms in Mol2 file %s (%i) does not\n",
              tfile->filename,frameAtom);
      mprinterr("       match number in associated parmtop (%i)!\n",natom);
      return -1;
    }
    if (frameAtom!=mol2atom) {
      mprintf("Warning: # atoms in Mol2 file (%s) frame %i (%i) not equal\n",
                tfile->filename,Frames,frameAtom);
      mprintf("         to # atoms int first frame (%i).\n",mol2atom);
      mprintf("         Only using frame 1.\n");
      Frames=1;
      break;
    }
    Frames++;
  }

  this->closeTraj();
  if (debug>0) mprintf("      Mol2 file %s has %i frames.\n",tfile->filename,Frames);
  seekable = false;
  hasTemperature=false;
  hasBox=false;

  return Frames;
}

/* Mol2File::readFrame()
 */
int Mol2File::readFrame(int set,double *X, double *V,double *box, double *T) {
  char buffer[MOL2BUFFERSIZE];
  int atom, atom3;

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(this->tfile, ATOM)) return 1;

  atom3=0;
  for (atom = 0; atom < mol2atom; atom++) {
    if (tfile->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    sscanf(buffer,"%*i %*s %lf %lf %lf",X+atom3, X+atom3+1, X+atom3+2);
    //F->printAtomCoord(atom);
    atom3+=3;
  }

  return 0;
}

/* Mol2File::SetWriteMode()
 * Set write mode to SINGLE, MOL, or MULTI
 */
void Mol2File::SetWriteMode(MOL2WRITEMODE modeIn) {
  mol2WriteMode = modeIn;
  //mprintf("MOL2 WRITE MODE SET TO %i\n",(int)mol2WriteMode);
}

/* Mol2File::NumFramesToWrite()
 * Modify write mode based on expected number of frames being written.
 */
void Mol2File::NumFramesToWrite(int parmFrames) {
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use @<TRIPOS>MOLECULE keyword to separate frames.
  if (mol2WriteMode==SINGLE && parmFrames>1) mol2WriteMode=MOL;
}

/* Mol2File::SetParmInfo()
 * In order to write mol2 files, need some parm info like atom names etc.
 * This should be called before setupWrite.
 */
void Mol2File::SetParmInfo(int nres, int nbondsh, int nbonds, NAME *names,
                           NAME *resnames, NAME *types, int *resnums,
                           int *bonds, int *bondsh, double *charge) {
  //mprintf("SETTING PARM INFO FOR Mol2 FILE\n"); // DEBUG
  trajnres = nres;
  trajnbondsh = nbondsh;
  trajnbonds = nbonds;
  trajAtomNames = names;
  trajTypes = types;
  // If types are not set just use atom names
  if (trajTypes == NULL) trajTypes = trajAtomNames;
  trajResNames = resnames;
  trajResNums = resnums;
  trajCharges = charge;
  trajBonds = bonds;
  trajBondsh = bondsh;
  // DEBUG
  //mprintf("Atomname 0 %s\n",trajAtomNames[0]);
  //mprintf("Resname 0  %s\n",trajResNames[0]);
  //mprintf("AtPerMol 0 %i\n",trajAtomsPerMol[0]);
  //mprintf("ResNums 0  %i\n",trajResNums[0]);
  //mprintf("Charge 0   %lf\n",trajCharges[0]);
}

/* Mol2File::setupWrite()
 * SetParmInfo should be called before calling this routine.
 */
int Mol2File::setupWrite(int natom) {
  mol2atom = natom;
  // Check number of bonds
  mol2bonds = trajnbondsh + trajnbonds;
  if (mol2bonds < 0) mol2bonds=0;
  if (mol2bonds == 0)
    mprintf("Warning: %s: topology does not contain bond information.\n",tfile->filename);
  // If more than 99999 atoms this may fail
  if (mol2atom>99999) {
    mprintf("Warning: %s: Large # of atoms (%i > 99999) for mol2 format.\n",
            tfile->filename,mol2atom);
    mprintf("         File may not write correctly.\n");
  }
  // Check that SetParmInfo has indeed been called
  if (trajAtomNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Atom names are NULL.\n",tfile->filename);
    return 1;
  }
  if (trajTypes==NULL) {
    mprinterr("Error: setupWrite [%s]: Atom types are NULL.\n",tfile->filename);
    return 1;
  }
  if (trajResNames==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue names are NULL.\n",tfile->filename);
    return 1;
  }
  if (trajResNums==NULL) {
    mprinterr("Error: setupWrite [%s]: Residue #s are NULL.\n",tfile->filename);
    return 1;
  }
  //if (trajCharges==NULL) {
  //  mprinterr("Error: setupWrite [%s]: Charges are NULL.\n",tfile->filename);
  //  return 1;
  //}
  if (mol2bonds > 0) {
    if (trajBonds==NULL) {
      mprinterr("Error: setupWrite [%s]: Bonds are NULL.\n",tfile->filename);
      return 1;
    }
    if (trajBondsh==NULL) {
      mprinterr("Error: setupWrite [%s]: BondsH are NULL.\n",tfile->filename);
      return 1;
    }
  }
  return 0;
}

/* Mol2File::writeFrame()
 */
int Mol2File::writeFrame(int set, double *X, double *V,double *box, double T) {
  char buffer[1024];
  int atom, atom3, res;
  double Charge;

  //mprintf("DEBUG: Calling Mol2File::writeFrame for set %i\n",set);
  if (mol2WriteMode==MULTI) {
    NumberFilename(buffer,tfile->filename,set+OUTPUTFRAMESHIFT);
    if (tfile->IO->Open(buffer,"wb")) return 1;
  }
  //@<TRIPOS>MOLECULE section
  tfile->IO->Printf("@<TRIPOS>MOLECULE\n");
  // mol_name
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  // mol_type
  // charge_type
  // [status_bits
  // [mol_comment]]
  tfile->IO->Printf("%s\n",title);
  tfile->IO->Printf("%5i %5i %5i %5i %5i\n",mol2atom,mol2bonds,1,0,0);
  tfile->IO->Printf("SMALL\n"); // May change this later
  if (trajCharges!=NULL)
    tfile->IO->Printf("USER_CHARGES\n"); // May change this later
  else
    tfile->IO->Printf("NO_CHARGES\n");
  tfile->IO->Printf("\n\n");

  //@<TRIPOS>ATOM section
  Charge = 0.0;
  tfile->IO->Printf("@<TRIPOS>ATOM\n");
  atom3=0;
  res = 0;
  for (atom=0; atom < mol2atom; atom++) {
    // figure out the residue number
    if ( atom==trajResNums[res+1] ) res++;
    if (trajCharges!=NULL) Charge = trajCharges[atom]; 
    tfile->IO->Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
                     atom+1, trajAtomNames[atom], X[atom3], X[atom3+1], X[atom3+2],
                     trajTypes[atom], res+1, trajResNames[res], Charge);
    atom3+=3;
  }

  //@<TRIPOS>BOND section
  if (mol2bonds > 0) {
    // Atom #s in the bonds and bondh array are * 3
    tfile->IO->Printf("@<TRIPOS>BOND\n");
    atom3=0;
    for (atom=0; atom < trajnbonds; atom++) {
      tfile->IO->Printf("%5d %5d %5d 1\n",atom+1,(trajBonds[atom3]/3)+1,(trajBonds[atom3+1]/3)+1);
      atom3+=3;
    }
    atom3=0;
    res = atom; // Where bonds left off
    for (atom=0; atom < trajnbondsh; atom++) {
      tfile->IO->Printf("%5d %5d %5d 1\n",res+1,(trajBondsh[atom3]/3)+1,(trajBondsh[atom3+1]/3)+1);
      atom3+=3;
      res++;
    }
  }

  //@<TRIPOS>SUBSTRUCTURE section
  tfile->IO->Printf("@<TRIPOS>SUBSTRUCTURE\n");
  for (res = 0; res < trajnres; res++) {
    tfile->IO->Printf("%7d %4s %14d ****               0 ****  **** \n",
                      res+1, trajResNames[res],trajResNums[res]+1);
  }

  // If writing 1 pdb per frame, close output file
  if (mol2WriteMode==MULTI)
    tfile->IO->Close();

  return 0;
}
 
/* Mol2File::info()
 */
void Mol2File::info() {
  mprintf("is a Tripos Mol2 file");
  if (tfile->access==WRITE) {
    if (mol2WriteMode==MULTI)
      mprintf(" (1 file per frame)");
    else if (mol2WriteMode==MOL)
      mprintf(" (1 MOLECULE per frame)");
  }
}
