// PDBfile
#include <cstdlib>
#include <cstring>
#include "PDBfile.h"
#include "PDBfileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
PDBfile::PDBfile() {
  pdbAtom=0;
}

// DESTRUCTOR
PDBfile::~PDBfile() { 
  //fprintf(stderr,"PDBfile Destructor.\n");
}
//------------------------------------------------------------------------
/* 
 * PDBfile::close()
 */
void PDBfile::close() {
  File->CloseFile();
}

/* 
 * PDBfile::open()
 */
int PDBfile::open() {

  if (File->OpenFile()) return 1;

  return 0;
}


/* 
 * PDBfile::SetupRead()
 * Scan PDB file to determine number of frames (models).
 */
int PDBfile::SetupRead() {
  int atom, scanPDB;

  if ( this->open() ) return 1;

  // Two strats - check for MODEL keywords or see how many times natom ATOMs can be read
  Frames=0;
  scanPDB=1;
  while (scanPDB) {
    atom=0;
    while ( atom < P->natom ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( File->IO->Gets(buffer,256) ) { scanPDB=0; break; }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      // Skip non-ATOM records
      if (strncmp(buffer,"ATOM",4)!=0 &&
          strncmp(buffer,"HETATM",6)!=0 ) continue;
      atom++;
    }
    if (scanPDB) Frames++;
  }
  this->close();

  if (Frames<1) {
    mprintf("ERROR: PDBfile::SetupRead(): No frames read. atom=%i expected %i.\n",
            atom,P->natom);
    return 1;
  }
  if (debug>0) mprintf("PDBfile::SetupRead(): %s %i atoms %i frames.\n",trajfilename,
                       atom,Frames);
  stop=Frames;
  pdbAtom = P->natom;
  return 0;
}

/* 
 * PDBfile::getFrame()
 * Read frame (model) from PDB file. Use pdbAtom from SetupRead instead of
 * P->natom in case of stripped prmtop.
 */
int PDBfile::getFrame(int set) {
  int atom, atom3;

  atom=0;
  atom3=0;
  while (atom < pdbAtom) {
    if ( File->IO->Gets(buffer,256) ) return 1;
    // Skip non-ATOM records
    if (strncmp(buffer,"ATOM",4)!=0 &&
        strncmp(buffer,"HETATM",6)!=0 ) continue;  
    // Read current PDB record XYZ into Frame
    //atom3 = atom * 3;
    pdb_xyz(buffer,F->X+atom3);
    atom++; 
    atom3+=3;
  }

  return 0;
}

/*
 * PDBfile::SetupWrite
 */ 
int PDBfile::SetupWrite( ) {
  return 0;
}

/*
 * PDBfile::writeFrame()
 * Write the frame (model) to PDB file.
 * NOTE: Eventually give option to write individual files or models.
 */
int PDBfile::writeFrame(int set) {
  int i,i3,res;
  float Occ, B;

  res=0; Occ=0.0; B=0.0;
  // Use F->natom instead of P->natom in case of stripped coordinates?
  i3=0;
  for (i=0; i<F->natom; i++) {
    // figure out the residue number
    if ( i==P->resnums[res+1] ) res++; 
    pdb_write_ATOM(buffer,"ATOM",i+1,P->names[i],P->ResidueName(res),'X',res+1,
                   F->X[i3],F->X[i3+1],F->X[i3+2],Occ,B,(char*)"\0");
    File->IO->Write(buffer,sizeof(char),strlen(buffer)); 
    i3+=3;
  }
  return 0;
}

/*
 * Info()
 */
void PDBfile::Info() {
  mprintf("  File (%s) is a PDB file", File->filename);
/*    if (p->option2 == 1) 
      printfone(" with no atom wrapping");
    if (p->option1 == 1)
      printfone(": AMBER charges and radii in prmtop to occupancy and temp factor columns");
    else if (p->option1 == 2)
      printfone(": AMBER charges and PARSE radii to occupancy and temp factor columns");
    else if (p->option1 == 3)
      printfone(": AMBER charges and vdw radii (r*) to occupancy and temp factor columns");*/
}
