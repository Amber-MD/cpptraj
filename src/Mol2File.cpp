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
  if (File->OpenFile()) return 1;

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
  mprintf("    Mol2 Title: [%s]\n",title);
  //   Scan # atoms
  if ( File->IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  sscanf(buffer,"%i",&mol2atom);
  mprintf("    Mol2 #atoms: %i\n",mol2atom);
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
    F->printAtomCoord(atom);
    atom3+=3;
  }

  return 0;
}

/*
 * Info()
 */
void Mol2File::Info() {
  mprintf("  File (%s) is a Tripos Mol2 file", File->filename);
}
