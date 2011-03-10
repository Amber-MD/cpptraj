// PDBfile
#include <cstdio> //sprintf
#include <cstdlib>
#include <cstring>
#include "PDBfile.h"
#include "PDBfileRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
PDBfile::PDBfile() {
  pdbAtom=0;
  writeMode=0;
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
  // Only close if not writing 1 pdb per frame
  if (writeMode!=2) {
    File->IO->Printf("%-6s\n","END");
    File->CloseFile();
  }
}

/* 
 * PDBfile::open()
 */
int PDBfile::open() {
  int err;

  err = 0; 
  switch (File->access) {
    case READ  : err = File->OpenFile(); break;
    case WRITE :
      // If writing 1 pdb per frame do not open here
      if (writeMode!=2) err = File->OpenFile(); 
      break;
    case APPEND:
      mprintf("Error: Append not supported for PDB files.\n");
      err=1;
      break;
  }
  
  return err;
}


/* 
 * PDBfile::SetupRead()
 * Scan PDB file to determine number of frames (models). The first frame will
 * also be checked to ensure that the atom names match those in the parm file.
 */
int PDBfile::SetupRead() {
  int atom, scanPDB;
  AmberParm::NAME namebuffer;

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
      // If still on frame 1, check to see if atom names match those in parm
      if (Frames==0) {
        pdb_name(buffer,namebuffer);
        if ( strcmp(P->names[atom], namebuffer)!=0 ) {
          mprintf("Warning: %s: Atom %i name [%s] does not match parm name [%s]\n",
                  trajfilename,atom,namebuffer,P->names[atom]);
        }
      }
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
 * PDBFile::WriteArgs()
 * Process arguments related to PDB write.
 */
int PDBfile::WriteArgs(ArgList *A) {
  // model: Multiple frames use model keyword (default if #frames<=1)
  if (A->hasKey("model")) writeMode=1; 
  // multi: Each frame is written to a different file
  if (A->hasKey("multi")) writeMode=2; 
  return 0;
}

/*
 * PDBfile::SetupWrite
 */ 
int PDBfile::SetupWrite( ) {
  // If writing more than 1 frame and not writing 1 pdb per frame, 
  // use MODEL keyword to separate frames.
  if (writeMode==0 && P->parmFrames>1) writeMode=1;
  return 0;
}

/*
 * PDBfile::writeFrame()
 * Write the frame (model) to PDB file.
 * NOTE: Eventually give option to write individual files or models.
 */
int PDBfile::writeFrame(int set) {
  int i,i3,res,atom,mol,lastAtomInMol;
  float Occ, B;

  // If writing 1 pdb per frame set up output filename and open
  if (writeMode==2) {
    sprintf(buffer,"%s.%i",File->filename,set + OUTPUTFRAMESHIFT);
    if (File->IO->Open(buffer,"wb")) return 1;
  // If specified, write MODEL keyword
  } else if (writeMode==1) {
    // 1-6 MODEL, 11-14 model serial #
    // Since num frames could be large, do not format the integer with width - OK?
    File->IO->Printf("MODEL     %i\n",set);
  }

  res=0; Occ=0.0; B=0.0;
  // Use F->natom instead of P->natom in case of stripped coordinates?
  i3=0;
  atom=1; // Actual PDB atom number
  mol=0;
  if (P->atomsPerMol!=NULL)
    lastAtomInMol=P->atomsPerMol[0];
  for (i=0; i<F->natom; i++) {
    // If this atom belongs to a new molecule print a TER card
    /*if (P->atomsPerMol!=NULL) {
      if (i == lastAtomInMol) {
        pdb_write_ATOM(buffer,PDBTER,atom,(char*)"",P->ResidueName(res),'X',res+1,
                       0,0,0,0,0,(char*)"\0");
        File->IO->Write(buffer,sizeof(char),strlen(buffer));
        atom++;
        mol++;
        if (mol<P->molecules) lastAtomInMol += P->atomsPerMol[mol];
      }
    }*/
    // figure out the residue number
    if ( i==P->resnums[res+1] ) res++; 
    pdb_write_ATOM(buffer,PDBATOM,atom,P->names[i],P->ResidueName(res),'X',res+1,
                   F->X[i3],F->X[i3+1],F->X[i3+2],Occ,B,(char*)"\0");
    File->IO->Write(buffer,sizeof(char),strlen(buffer)); 
    i3+=3;
    atom++;
  }

  // If writing 1 pdb per frame, close output file
  if (writeMode==2) {
    File->IO->Printf("%-6s\n","END");
    File->IO->Close();
  // If MODEL keyword was written, write corresponding ENDMDL record
  } else if (writeMode==1) {
    File->IO->Printf("ENDMDL\n");
  }
  return 0;
}

/*
 * Info()
 */
void PDBfile::Info() {
  mprintf("is a PDB file");
  if (File->access==WRITE) {
    if (writeMode==2)
      mprintf(" (1 file per frame)");
    else if (writeMode==1)
      mprintf(" (1 MODEL per frame)");
  }
   
/*    if (p->option2 == 1) 
      printfone(" with no atom wrapping");
    if (p->option1 == 1)
      printfone(": AMBER charges and radii in prmtop to occupancy and temp factor columns");
    else if (p->option1 == 2)
      printfone(": AMBER charges and PARSE radii to occupancy and temp factor columns");
    else if (p->option1 == 3)
      printfone(": AMBER charges and vdw radii (r*) to occupancy and temp factor columns");*/
}
