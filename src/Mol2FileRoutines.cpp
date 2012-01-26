#include <cstring>
#include <cstdio>
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
/*! \file Mol2FileRoutines.cpp
    \brief Collection of routines used to access mol2 files.
 */

/// Tripos Tags - must be in same order as enum type TRIPOSTAG
const char TRIPOSTAGTEXT[NUMTRIPOSTAGS][30]={
  "@<TRIPOS>MOLECULE\0",
  "@<TRIPOS>ATOM\0",
  "@<TRIPOS>BOND\0",
  "@<TRIPOS>SUBSTRUCTURE\0"
};

// Mol2ScanTo()
/** Scan to the specified TRIPOS section of file
  */
int Mol2ScanTo( CpptrajFile *File, TRIPOSTAG tag ) {
  int tagSize;
  char buffer[MOL2BUFFERSIZE];

  tagSize = strlen(TRIPOSTAGTEXT[tag]);
  while ( File->IO->Gets(buffer,MOL2BUFFERSIZE)==0 ) {
    //mprintf("DEBUG: Line [%s]\n",buffer);
    //mprintf("DEBUG: Targ [%s]\n",TRIPOSTAGTEXT[tag]); 
    if (strncmp(buffer,TRIPOSTAGTEXT[tag],tagSize)==0) return 0;
  }

  // Suppress this warning so routine can be used to scan # frames
  //mprintf("Warning: Mol2File::ScanTo(): Could not find tag %s\n",TRIPOSTAGTEXT[tag]);

  return 1;
}

// Mol2AtomName()
/** Given a Mol2 ATOM line, return atom name. Trim to 4 chars to be consistent 
  * with the rest of Amber. 
  */
int Mol2AtomName(char *buffer, NAME name) {
  char tmp[10];
  if (buffer==NULL || name==NULL) return 1;
  sscanf(buffer,"%*i %s",tmp);
  name[0]=tmp[0];
  name[1]=tmp[1];
  name[2]=tmp[2];
  name[3]=tmp[3];
  PadWithSpaces(name); 
  name[4]='\0';
  // Replace asterisks with prime to prevent atom mask problems
  ReplaceAsterisk(name);
  //mprintf("DEBUG: MOL2: name [%s]\n",name);
  return 0;
}

// Mol2XYZ()
/** Given a Mol2 ATOM line, get the X Y and Z coords.
  */
int Mol2XYZ(char *buffer, double *X) {
  if (buffer==NULL || X==NULL) return 1;
  sscanf(buffer,"%*i %*s %lf %lf %lf",X, X+1, X+2);
  return 0;
}

// Mol2AtomType
/** Given a Mol2 ATOM line, return atom type. Try to convert Sybyl atom type
  * to amber type.
  * Sybyl atom types seem to have at most 5 chars - increasing size of NAME
  * huld be ok.
  */
int Mol2AtomType(char *buffer, NAME type) {
  char tmp[10];
  if (buffer==NULL || type==NULL) return 1;
  sscanf(buffer,"%*i %*s %*f %*f %*f %s",tmp);
/*
  if      (strcmp(tmp,"C.3")==0)   strcpy(type,"CT  ");
  else if (strcmp(tmp,"C.2")==0)   strcpy(type,"C   ");
  else if (strcmp(tmp,"C.1")==0)   strcpy(type,"CZ  ");
  else if (strcmp(tmp,"C.ar")==0)  strcpy(type,"CA  ");
  else if (strcmp(tmp,"C.cat")==0) strcpy(type,"C+  ");
  else if (strcmp(tmp,"N.3")==0)   strcpy(type,"NT  ");
  else if (strcmp(tmp,"N.2")==0)   strcpy(type,"N2  ");
  else if (strcmp(tmp,"N.1")==0)   strcpy(type,"NY  ");
  else if (strcmp(tmp,"N.ar")==0)  strcpy(type,"NC  ");
  else if (strcmp(tmp,"N.am")==0)  strcpy(type,"N   ");
  else if (strcmp(tmp,"N.pl3")==0) strcpy(type,"N   "); // NOTE: Nitro, not in ff99
  else if (strcmp(tmp,"N.4")==0)   strcpy(type,"N3  ");
  else if (strcmp(tmp,"O.3")==0)   strcpy(type,"OS  "); // ??
  else if (strcmp(tmp,"O.2")==0)   strcpy(type,"O   "); 
  else if (strcmp(tmp,"O.co2")==0) strcpy(type,"O2  "); 
  else if (strcmp(tmp,"O.spc")==0) strcpy(type,"OW  "); 
  else if (strcmp(tmp,"O.t3p")==0) strcpy(type,"OW  "); 
  else if (strcmp(tmp,"S.3")==0)   strcpy(type,"SH  "); 
  else if (strcmp(tmp,"S.2")==0)   strcpy(type,"S   "); 
  else if (strcmp(tmp,"S.O")==0)   strcpy(type,"S   "); 
  else if (strcmp(tmp,"S.O2")==0)  strcpy(type,"S   "); 
  else if (strcmp(tmp,"P.3")==0)   strcpy(type,"P   "); 
  else if (strcmp(tmp,"H.spc")==0) strcpy(type,"HW  "); 
  else if (strcmp(tmp,"H.t3p")==0) strcpy(type,"HW  ");
  else {*/
    // Copy over the first 5 chars
    strncpy(type,tmp,5);
    type[5]='\0';
  //}
  return 0;
}

// Mol2ResNumName()
/** Get residue number and name, trim to 4 chars
  */
int Mol2ResNumName(char *buffer, int *resnum, NAME resname) {
  char tmp[10];
  if (buffer==NULL || resnum==NULL || resname==NULL) return 1;
  sscanf(buffer,"%*i %*s %*f %*f %*f %*s %i %s",resnum,tmp);
  resname[0]=tmp[0]; 
  resname[1]=tmp[1]; 
  resname[2]=tmp[2]; 
  resname[3]=tmp[3];
  PadWithSpaces(resname);
  resname[4]='\0';
  // Replace asterisks with prime to prevent atom mask problems
  ReplaceAsterisk(resname);
  return 0;
}

// Mol2Charge()
/** Get charge */
double Mol2Charge(char *buffer) {
  double q;
  sscanf(buffer,"%*i %*s %*f %*f %*f %*s %*i %*s %lf",&q);
  return q;
} 

