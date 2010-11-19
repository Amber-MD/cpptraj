// PDBfile
#include <cstdlib>
#include <cstring>
#include "PDBfile.h"
#include "PDBfileRoutines.h"

// CONSTRUCTOR
PDBfile::PDBfile() { }

// DESTRUCTOR
PDBfile::~PDBfile() { 
  //fprintf(stderr,"PDBfile Destructor.\n");
}
//------------------------------------------------------------------------
/* PDBfile::close()
 */
void PDBfile::close() {
  File->CloseFile();
}

/* PDBfile::open()
 */
int PDBfile::open() {

  if (File->OpenFile()) return 1;

  return 0;
}


/* PDBfile::SetupRead()
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
    fprintf(stdout,"ERROR: PDBfile::SetupRead(): No frames read. atom=%i expected %i.\n",
            atom,P->natom);
    return 1;
  }
  //seekable=1; // Set seekable since only 1 frame read (i.e. we know size)
  //Frames=1;
  if (debug>0) fprintf(stdout,"PDBfile::SetupRead(): %s %i atoms %i frames.\n",trajfilename,
                       atom,Frames);
  stop=Frames;
  return 0;
}

/* PDBfile::getFrame()
 */
int PDBfile::getFrame(int set) {
  int atom, atom3;

  atom=0;

  while (atom < P->natom) {
    if ( File->IO->Gets(buffer,256) ) return 1;
    // Skip non-ATOM records
    if (strncmp(buffer,"ATOM",4)!=0 &&
        strncmp(buffer,"HETATM",6)!=0 ) continue;  
    // Read current PDB record XYZ into Frame
    atom3 = atom * 3;
    pdb_xyz(buffer,F->X+atom3); 
    atom++;
  }

  return 0;
}

// Set up trajectory for either write or append
int PDBfile::SetupWrite( ) {
  return 0;
}

// Write a frame
int PDBfile::writeFrame(int set) {
  int i,i3,res, resOut;
  double Occ, B;

  res=0; resOut=0; Occ=0.0; B=0.0;
  for (i=0; i<P->natom; i++) {
    i3 = i * 3;
    // figure out the residue number
    if ( (i+1)==P->resnums[res+1] ) {res++; resOut++;}
    if (resOut>9998) resOut=0;
    File->IO->Printf("%-6s%5i %-4s%-4s %c%4i    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf%14s\n",
            "ATOM",i+1,P->names[i],P->resnames[res],'X',resOut+1,
            F->X[i3],F->X[i3+1],F->X[i3+2],Occ,B,"");
  }
  return 0;
}

/*
 * Info()
 */
void PDBfile::Info() {
  fprintf(stdout,"  File (%s) is a PDB file", File->filename);
/*    if (p->option2 == 1) 
      printfone(" with no atom wrapping");
    if (p->option1 == 1)
      printfone(": AMBER charges and radii in prmtop to occupancy and temp factor columns");
    else if (p->option1 == 2)
      printfone(": AMBER charges and PARSE radii to occupancy and temp factor columns");
    else if (p->option1 == 3)
      printfone(": AMBER charges and vdw radii (r*) to occupancy and temp factor columns");*/
}
