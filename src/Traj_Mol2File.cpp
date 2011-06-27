// Mol2File
#include <cstdio> // sscanf
#include "Traj_Mol2File.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Mol2File::Mol2File() {
  mol2atom=0;
  mol2bonds=0;
  Types=NULL;
  mol2WriteMode=SINGLE;
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
int Mol2File::readFrame(int set,double *X, double *box, double *T) {
  char buffer[MOL2BUFFERSIZE];
  int atom, atom3;

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(this->File, ATOM)) return 1;

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

/* Mol2File::setupWrite()
 * SetParmInfo should be called before calling this routine.
 */
int Mol2File::setupWrite() {
  // Check number of bonds
  mol2bonds = P->NbondsWithH() + P->NbondsWithoutH();
  if (mol2bonds < 0) mol2bonds=0;
  if (mol2bonds == 0)
    mprintf("Warning: Mol2File::SetupWrite: %s does not contain bond information.\n",
            P->parmName);
  // If more than 99999 atoms this will fail
  if (P->natom>99999) {
    mprintf("Warning: Mol2File::SetupWrite: %s: too many atoms for mol2 format.\n",
            P->parmName);
    mprintf("         File %s may not write correctly.\n",trajfilename);
  }
  // If types are not set just use atom names
  Types = P->types;
  if (P->types==NULL)
    Types = P->names;
  // If writing more than 1 frame and not writing 1 mol2 per frame,
  // use MOLECULE to separate frames.
  if (writeMode==0 && P->parmFrames>1) writeMode=1;

  return 0;
}


/*
 * Mol2File::writeFrame()
 * NOTE: P->natom should equal F->natom!
 */
int Mol2File::writeFrame(int set) {
  char buffer[1024];
  int atom, atom3, res;
  double Charge;

  //mprintf("DEBUG: Calling Mol2File::writeFrame for set %i\n",set);
  if (writeMode==2) {
    sprintf(buffer,"%s.%i",tfile->filename,set+OUTPUTFRAMESHIFT);
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
  tfile->IO->Printf("%5i %5i %5i %5i %5i\n",F->natom,mol2bonds,1,0,0);
  tfile->IO->Printf("SMALL\n"); // May change this later
  if (P->charge!=NULL)
    tfile->IO->Printf("USER_CHARGES\n"); // May change this later
  else
    tfile->IO->Printf("NO_CHARGES\n");
  tfile->IO->Printf("\n\n");

  //@<TRIPOS>ATOM section
  Charge = 0.0;
  tfile->IO->Printf("@<TRIPOS>ATOM\n");
  atom3=0;
  for (atom=0; atom < F->natom; atom++) {
    res = P->atomToResidue(atom);
    if (P->charge!=NULL) Charge = P->charge[atom]; 
    tfile->IO->Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
                     atom+1, P->names[atom], F->X[atom3], F->X[atom3+1], F->X[atom3+2],
                     Types[atom], res+1, P->ResidueName(res), Charge);
    atom3+=3;
  }

  //@<TRIPOS>BOND section
  if (mol2bonds > 0) {
    // Atom #s in the bonds and bondh array are * 3
    tfile->IO->Printf("@<TRIPOS>BOND\n");
    atom3=0;
    for (atom=0; atom < P->NbondsWithoutH(); atom++) {
      tfile->IO->Printf("%5d %5d %5d 1\n",atom+1,(P->bonds[atom3]/3)+1,(P->bonds[atom3+1]/3)+1);
      atom3+=3;
    }
    atom3=0;
    res = atom; // Where bonds left off
    for (atom=0; atom < P->NbondsWithH(); atom++) {
      tfile->IO->Printf("%5d %5d %5d 1\n",res+1,(P->bondsh[atom3]/3)+1,(P->bondsh[atom3+1]/3)+1);
      atom3+=3;
      res++;
    }
  }

  //@<TRIPOS>SUBSTRUCTURE section
  tfile->IO->Printf("@<TRIPOS>SUBSTRUCTURE\n");
  for (res = 0; res < P->nres; res++) {
    tfile->IO->Printf("%7d %4s %14d ****               0 ****  **** \n",
                      res+1, P->ResidueName(res),P->resnums[res]+1);
  }

  // If writing 1 pdb per frame, close output file
  if (writeMode==2)
    tfile->IO->Close();

  currentFrame++;

  return 0;
}
 
/*
 * Info()
 */
void Mol2File::Info() {
  mprintf("is a Tripos Mol2 file");
  if (tfile->access==WRITE) {
    if (writeMode==2)
      mprintf(" (1 file per frame)");
    else if (writeMode==1)
      mprintf(" (1 MOLECULE per frame)");
  }
}
