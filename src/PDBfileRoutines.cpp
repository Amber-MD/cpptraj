/* PDBfileRoutines.cpp
 * A collection of functions that will parse a line from a PDB file
 * and return the requested information.
 */
#include <cstdlib>
#include <cstring>
#include <cstdio> //sprintf
#include <cctype> //isdigit
#include "PDBfileRoutines.h"

// PDB record types
const char PDB_RECNAME[3][7]={"ATOM", "HETATM", "TER"};

// ====================== PRIVATE ROUTINES ===================================
void TrimName(char *);
void WrapName(char *);
void ReplaceAsterisk(char *);
/*
 * TrimName()
 * Given a string of length 5 (4 chars + 1 NULL) remove leading whitespace
 */
void TrimName(char *NameIn) {
  if        (NameIn[0]!=' ') { // No leading whitespace
    return;
  } else if (NameIn[1]!=' ') { // [_XXX]
    NameIn[0]=NameIn[1];
    NameIn[1]=NameIn[2];
    NameIn[2]=NameIn[3];
    NameIn[3]=' ';
  } else if (NameIn[2]!=' ') { // [__XX]
    NameIn[0]=NameIn[2];
    NameIn[1]=NameIn[3];
    NameIn[2]=' ';
    NameIn[3]=' ';
  } else if (NameIn[3]!=' ') { // [___X]
    NameIn[0]=NameIn[3];
    NameIn[1]=' ';
    NameIn[2]=' ';
    NameIn[3]=' ';
  }
  // Otherwise Res is Blank, no trim needed
  return;
}

/*
 * WrapName()
 * Given a string of length 5 (4 chars + 1 NULL), move leading numbers to
 * the back of the string until first char is alphabetic.
 */
void WrapName(char *NameIn) {
  int i;
  int blank=-1;
  int numalpha=0;
  char digit;
  //fprintf(stderr,"NameIn = [%s]\n",NameIn);
  // First trim leading whitespace if necessary
  // NOTE: Should already be done
  //TrimName(NameIn);
  // Check for last blank char and alpha chars
  for (i=0; i<4; i++) {
    if (!isdigit(NameIn[i])) {
      // If the first char is already alpha exit
      if (i==0) return;
      numalpha++;
    }
    if (blank<0 && NameIn[i]==' ') blank=i;
  }
  if (blank==-1) blank=4;
  // No alpha chars, return
  if (numalpha==0) return;
  // No chars (____), just return
  if (blank==0) return;
  // 1 char (X___), also return 
  if (blank==1) return;
  // Otherwise 2 char (XX__), 3 char (XXX_), or 4 char (XXXX)
  while (isdigit(NameIn[0])) {
    digit=NameIn[0];
    for (i=1; i<blank; i++) NameIn[i-1]=NameIn[i];
    NameIn[blank-1]=digit;
  }
  //fprintf(stderr,"NameOut= [%s]\n",NameIn);
}

/*
 * ReplaceAsterisk()
 * Given a string of length 5 (4chars + 1 NULL) change any asterisk (*) to
 * prime ('). In cpptraj asterisks are considered reserved characters for
 * atom masks.
 */
void ReplaceAsterisk(char *NameIn) {
  if (NameIn[0]=='*') NameIn[0]='\'';
  if (NameIn[1]=='*') NameIn[1]='\'';
  if (NameIn[2]=='*') NameIn[2]='\'';
  if (NameIn[3]=='*') NameIn[3]='\'';
  return;
}
// ===========================================================================

/*
 * isPDBkeyword()
 * Return true if the first 6 chars of buffer match a PDB keyword
 */
bool isPDBkeyword(char *buffer) {
  if (strncmp(buffer,"HEADER",6)==0) return true;
  if (strncmp(buffer,"TITLE ",6)==0) return true;
  if (strncmp(buffer,"COMPND",6)==0) return true;
  if (strncmp(buffer,"AUTHOR",6)==0) return true;
  if (strncmp(buffer,"ATOM  ",6)==0) return true;
  if (strncmp(buffer,"HETATM",6)==0) return true;
  if (strncmp(buffer,"CRYST1",6)==0) return true;
  if (strncmp(buffer,"REMARK",6)==0) return true;
  if (strncmp(buffer,"MODEL ",6)==0) return true;
  if (strncmp(buffer,"JRNL  ",6)==0) return true;
  if (strncmp(buffer,"SEQRES",6)==0) return true;
  return false;
}

// PDB Record Title
char *pdb_title(char *buffer) {
  char *title;

  title=(char*) malloc(7*sizeof(char));
  title[0]=buffer[0];
  title[1]=buffer[1];
  title[2]=buffer[2];
  title[3]=buffer[3];
  title[4]=buffer[4];
  title[5]=buffer[5];
  title[6]='\0';

  return title;
}

// PDB Record Atom Number
int pdb_atom(char *buffer) {
  char temp[6];

  temp[0]=buffer[6];
  temp[1]=buffer[7];
  temp[2]=buffer[8];
  temp[3]=buffer[9];
  temp[4]=buffer[10];
  temp[5]='\0';
  return atoi(temp);
}

// PDB Record Atom Name
int pdb_name(char *buffer, char *name) {
  name[0]=buffer[12];
  name[1]=buffer[13];
  name[2]=buffer[14];
  name[3]=buffer[15];
  name[4]='\0';
  // Trim Leading whitespace
  TrimName(name);
  // Wrap name if it starts with a digit
  WrapName(name);
  // Replace asterisks with prime to prevent atom mask problems
  ReplaceAsterisk(name);
  return 0;
}

// PDB Record Residue Name
int pdb_resname(char *buffer, char *resname) {
  resname[0]=buffer[16]; // Alternate location indicator
  resname[1]=buffer[17];
  resname[2]=buffer[18];
  resname[3]=buffer[19];
  resname[4]='\0';
  // Trim Leading whitespace 
  TrimName(resname);
  // Replace asterisks with prime to prevent atom mask problems
  ReplaceAsterisk(resname);
  return 0;
}

// PDB Record Chain ID
char pdb_chain(char *buffer) {
  return buffer[21];
}

// PDB Record Residue Number
int pdb_resnum(char *buffer) {
  char temp[6];

  temp[0]=buffer[22];
  temp[1]=buffer[23];
  temp[2]=buffer[24];
  temp[3]=buffer[25];
  temp[4]=buffer[26]; // Code for insertion of residues
  temp[5]='\0';

  return atoi(temp);
} 

// PDB X Y and Z
// Memory for double array should be allocated
int pdb_xyz(char *buffer, double *X) {
  char temp[9];

  temp[0]=buffer[30];
  temp[1]=buffer[31];
  temp[2]=buffer[32];
  temp[3]=buffer[33];
  temp[4]=buffer[34];
  temp[5]=buffer[35];
  temp[6]=buffer[36];
  temp[7]=buffer[37];
  temp[8]='\0';
  X[0]=atof(temp);

  temp[0]=buffer[38];
  temp[1]=buffer[39];
  temp[2]=buffer[40];
  temp[3]=buffer[41];
  temp[4]=buffer[42];
  temp[5]=buffer[43];
  temp[6]=buffer[44];
  temp[7]=buffer[45];
  temp[8]='\0';
  X[1]=atof(temp);

  temp[0]=buffer[46];
  temp[1]=buffer[47];
  temp[2]=buffer[48];
  temp[3]=buffer[49];
  temp[4]=buffer[50];
  temp[5]=buffer[51];
  temp[6]=buffer[52];
  temp[7]=buffer[53];
  temp[8]='\0';
  X[2]=atof(temp);

  return 0;
}

// PDB Record Occupancy 
double pdb_occ(char *buffer) {
  char temp[7];

  temp[0]=buffer[54];
  temp[1]=buffer[55];
  temp[2]=buffer[56];
  temp[3]=buffer[57];
  temp[4]=buffer[58];
  temp[5]=buffer[59];
  temp[6]='\0';
  return atof(temp);
}

// PDB Record B-factor
double pdb_Bfactor(char *buffer) {
  char temp[7];

  temp[0]=buffer[60];
  temp[1]=buffer[61];
  temp[2]=buffer[62];
  temp[3]=buffer[63];
  temp[4]=buffer[64];
  temp[5]=buffer[65];
  temp[6]='\0';
  return atof(temp);
}

// 10 chars between Bfactor and element
char *pdb_lastChar(char *buffer) {
  char *E;

  E=(char*) malloc(11*sizeof(char));
  E[0]=buffer[66];
  E[1]=buffer[67];
  E[2]=buffer[68];
  E[3]=buffer[69];
  E[4]=buffer[70];
  E[5]=buffer[71];
  E[6]=buffer[72];
  E[7]=buffer[73];
  E[8]=buffer[74];
  E[9]=buffer[75];
  E[10]='\0';
  return E;
}

// Element. If blank, try to guess from the name
// H C N O F Cl Br S P
char *pdb_elt(char *buffer) {
  char *E,*ptr;
  char name[5];

  E=(char*) malloc(3*sizeof(char));
  E[0]=buffer[76]; // Element
  E[1]=buffer[77]; // Element
  E[2]='\0';

  if (E[0]==' ' && E[1]==' ') {
    pdb_name(buffer,name);
    // position ptr at first non-space character
    ptr=name;
    while (*ptr==' ' && *ptr!='\0') ptr++;
    // if NULL something went wrong, abort
    if (*ptr=='\0') return E;
    E[0]=ptr[0];
    // If C, check for L or l for chlorine
    if (ptr[0]=='C') {
      if (ptr[1]=='L' || ptr[1]=='l') E[1]='l';
    }
    // If B, check for R or r for bromine
    if (ptr[0]=='B') {
      if (ptr[1]=='R' || ptr[1]=='r') E[1]='r';
    }
  }

  return E;
}

// Charge. 
char *pdb_charge(char *buffer) {
  char *E;

  E=(char*) malloc(3*sizeof(char));
  E[0]=buffer[78]; // Charge
  E[1]=buffer[79]; // Charge
  E[2]='\0';
  return E;
}

// Write out an ATOM or HETATM record
void pdb_write_ATOM(char *buffer, PDB_RECTYPE Record, int atom, char *name, 
                    char *resnameIn, char chain, int resnum, 
                    double X, double Y, double Z, float Occ, float B,
                    char *End, bool highPrecision) {
  char resName[5];
  char *ptr;

  resName[4]='\0';
  // Residue number in PDB format can only be 4 digits wide
  while (resnum>9999) resnum-=9999;
  // Residue names in PDB format are 3 chars long starting at column 18. 
  // However in Amber residues are 4 characters long, usually with a space
  // at the end. If this is the case remove the space so that the residue name
  // lines up properly.
  resName[0] = resnameIn[0];
  resName[1] = resnameIn[1];
  resName[2] = resnameIn[2];
  if (resnameIn[3]!=' ')
    resName[3] = resnameIn[3];
  else
    resName[3] = '\0';

  ptr=buffer;
  sprintf(ptr,"%-6s%5i %-4s%4s %c%4i",PDB_RECNAME[Record], atom, name, resName, chain, resnum);
  ptr+=26;
  if (Record == PDBTER) {
    sprintf(ptr,"\n");
    return;
  }
  if (highPrecision)
    sprintf(ptr,"    %8.3lf%8.3lf%8.3lf%8.4f%8.4f%10s\n", X, Y, Z, Occ, B, End);
  else
    sprintf(ptr,"    %8.3lf%8.3lf%8.3lf%6.2f%6.2f%14s\n", X, Y, Z, Occ, B, End);
}
                    
