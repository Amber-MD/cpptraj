#include <cstdio> // sscanf
#include "Mol2File.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// Mol2File

// CONSTRUCTOR
Mol2File::Mol2File() {
  mol2atom=0;
}

// DESTRUCTOR
Mol2File::~Mol2File() {
}

/*
 * Mol2File::open()
 */
int Mol2File::open() {

  switch (File->access) {
    case READ :
      if (File->OpenFile()) return 1;
      break;
    case APPEND :
      mprintf("Error: Append not supported for mol2 files.\n");
      return 1;
      break;
    case WRITE :
      // Set up title
      if (title==NULL)
        this->SetTitle((char*)"Cpptraj generated mol2 file.\0");
      break;
  }
      
  return 0;
}

/*
 * Mol2File::close() {
 */
void Mol2File::close() {
  File->CloseFile();
}

/*
 * Mol2File::SetupRead()
 */
int Mol2File::SetupRead() {
  char buffer[MOL2BUFFERSIZE];

  if (this->open()) return 1;

  // Get @<TRIPOS>MOLECULE information
  if (Mol2ScanTo(this->File, MOLECULE)) return 1;
  //   Scan title
  if ( File->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1; 
  this->SetTitle(buffer);
  if (debug>0) mprintf("    Mol2 Title: [%s]\n",title);
  //   Scan # atoms
  if ( File->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  sscanf(buffer,"%i",&mol2atom);
  if (debug>0) mprintf("    Mol2 #atoms: %i\n",mol2atom);
  if (mol2atom!=P->natom) {
    mprintf("Error: # atoms in mol2 file (%i) not equal to # atoms\n");
    mprintf("       in associated parm file %s (%i).\n",P->parmName,P->natom);
  }

  // Currently allow only 1 Frame
  Frames=1;
  stop=Frames;

  return 0;
}

/*
 * Mol2File::getFrame()
 */
int Mol2File::getFrame(int set) {
  char buffer[MOL2BUFFERSIZE];
  int atom, atom3;

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(this->File, ATOM)) return 1;

  atom3=0;
  for (atom = 0; atom < mol2atom; atom++) {
    if (File->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    sscanf(buffer,"%*i %*s %lf %lf %lf",F->X+atom3, F->X+atom3+1, F->X+atom3+2);
    //F->printAtomCoord(atom);
    atom3+=3;
  }

  return 0;
}

/*
 * Mol2File::SetupWrite()
 */
int Mol2File::SetupWrite() {
  return 0;
}


/*
 * Mol2File::writeFrame()
 */
int Mol2File::writeFrame(int set) {
  char buffer[1024];
  AmberParm::NAME *Types;
  int mol2bonds, atom, atom3, res;
  double Charge;

  //mprintf("DEBUG: Calling Mol2File::writeFrame for set %i\n",set);  
  sprintf(buffer,"%s.%i",File->filename,set+OUTPUTFRAMESHIFT);
  if (File->IO->Open(buffer,"wb")) return 1;

  //@<TRIPOS>MOLECULE section
  File->IO->Printf("@<TRIPOS>MOLECULE\n");
  // If more than 99999 atoms this will fail
  if (F->natom>99999) 
    mprintf("Warning: Too many atoms for mol2 format - file may not write correctly.\n");
  // mol_name
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  // mol_type
  // charge_type
  // [status_bits
  // [mol_comment]]
  mol2bonds = P->NbondsWithH() + P->NbondsWithoutH();
  if (mol2bonds < 0) mol2bonds=0;
  File->IO->Printf("%s\n",title);
  File->IO->Printf("%5i %5i %5i %5i %5i\n",F->natom,mol2bonds,1,0,0);
  File->IO->Printf("SMALL\n"); // May change this later
  if (P->charge!=NULL)
    File->IO->Printf("USER_CHARGES\n"); // May change this later
  else
    File->IO->Printf("NO_CHARGES\n");
  File->IO->Printf("\n\n");

  //@<TRIPOS>ATOM section
  // If types are not set just use atom names
  Types = P->types;
  if (P->types==NULL)
    Types = P->names;
  Charge = 0.0;
  File->IO->Printf("@<TRIPOS>ATOM\n");
  atom3=0;
  for (atom=0; atom < F->natom; atom++) {
    res = P->atomToResidue(atom);
    if (P->charge!=NULL) Charge = P->charge[atom]; 
    File->IO->Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-4s %6i %-6s %10.6lf\n",
                     atom+1, P->names[atom], F->X[atom3], F->X[atom3+1], F->X[atom3+2],
                     Types[atom], res+1, P->ResidueName(res), Charge);
    atom3+=3;
  }

  //@<TRIPOS>BOND section
  if (mol2bonds > 0) {
    // Atom #s in the bonds and bondh array are * 3
    File->IO->Printf("@<TRIPOS>BOND\n");
    atom3=0;
    for (atom=0; atom < P->NbondsWithoutH(); atom++) {
      File->IO->Printf("%5d %5d %5d 1\n",atom+1,(P->bonds[atom3]/3)+1,(P->bonds[atom3+1]/3)+1);
      atom3+=3;
    }
    atom3=0;
    res = atom; // Where bonds left off
    for (atom=0; atom < P->NbondsWithH(); atom++) {
      File->IO->Printf("%5d %5d %5d 1\n",res+1,(P->bondsh[atom3]/3)+1,(P->bondsh[atom3+1]/3)+1);
      atom3+=3;
      res++;
    }
  } else {
    mprintf("Warning: Mol2File::writeFrame: Parm does not contain bond information.\n");
  }

  //@<TRIPOS>SUBSTRUCTURE section
  File->IO->Printf("@<TRIPOS>SUBSTRUCTURE\n");
  for (res = 0; res < P->nres; res++) {
    File->IO->Printf("%7d %4s %14d ****               0 ****  **** \n",
                      res+1, P->ResidueName(res),P->resnums[res]+1);
  }

  File->IO->Close();

  return 0;
}
 
/*
 * Info()
 */
void Mol2File::Info() {
  mprintf("is a Tripos Mol2 file");
}
