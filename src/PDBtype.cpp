#include <cstdlib> // atoi, atof
#include <cstring> //strncmp, strlen
#include <cstdio>  //sprintf
#include "PDBtype.h"

PDBtype::PDBtype() {
  buffer_[0]='\0';
  XYZ_[0] = 0;
  XYZ_[1] = 0;
  XYZ_[2] = 0;
}

/// PDB record types
const char PDBtype::PDB_RECNAME[3][7]={"ATOM", "HETATM", "TER"};

bool PDBtype::PDB_GetNextRecord(FileIO *IO) {
  return ( IO->Gets(buffer_, BUF_SIZE_) == 0 );
}

// PDBtype::IsPDBkeyword()
/** \return true if the first 6 chars of buffer match a PDB keyword */
bool PDBtype::IsPDBkeyword() {
  if (strncmp(buffer_,"HEADER",6)==0) return true;
  if (strncmp(buffer_,"TITLE ",6)==0) return true;
  if (strncmp(buffer_,"COMPND",6)==0) return true;
  if (strncmp(buffer_,"AUTHOR",6)==0) return true;
  if (strncmp(buffer_,"ATOM  ",6)==0) return true;
  if (strncmp(buffer_,"HETATM",6)==0) return true;
  if (strncmp(buffer_,"CRYST1",6)==0) return true;
  if (strncmp(buffer_,"REMARK",6)==0) return true;
  if (strncmp(buffer_,"MODEL ",6)==0) return true;
  if (strncmp(buffer_,"JRNL  ",6)==0) return true;
  if (strncmp(buffer_,"SEQRES",6)==0) return true;
  return false;
}

// PDBtype::ID()
/** Check that the first two lines contain valid PDB records. */
bool PDBtype::ID(FileIO *IO) {
  if (!PDB_GetNextRecord(IO)) return false;
  if (!IsPDBkeyword()) return false;
  if (!PDB_GetNextRecord(IO)) return false;
  if (!IsPDBkeyword()) return false;
  return true;
}

// PDBtype::IsPDBatomKeyword
/** \return true if the first 6 chars match ATOM or HETATM */
bool PDBtype::IsPDBatomKeyword() {
  if (strncmp(buffer_,"ATOM  ",6)==0) return true;
  if (strncmp(buffer_,"HETATM",6)==0) return true;
  return false;
}

bool PDBtype::IsPDB_TER() {
  if (buffer_[0]=='T' &&
      buffer_[1]=='E' &&
      buffer_[2]=='R'   )
    return true;
  return false;
}

bool PDBtype::IsPDB_END() {
  if (buffer_[0]=='E' &&
      buffer_[1]=='N' &&
      buffer_[2]=='D'   )
    return true;
  return false;
}

// PDBtype::pdb_Atom()
Atom PDBtype::pdb_Atom() {
  // Atom number (6-11)
  // Atom name (12-16)
  char savechar = buffer_[16];
  buffer_[16] = '\0';
  NameType aname(buffer_+12);
  // Replace asterisks with single quotes
  aname.ReplaceAsterisk();
  buffer_[16] = savechar;
  // X coord (30-38)
  // Y coord (38-46)
  // Z coord (46-54)

  return Atom(aname);
}

// PDBtype::pdb_Residue()
Residue PDBtype::pdb_Residue() {
  // Res name (16-20)
  char savechar = buffer_[20];
  buffer_[20] = '\0';
  NameType resname(buffer_+16);
  // Replace asterisks with single quotes
  resname.ReplaceAsterisk();
  buffer_[20] = savechar;
  // Res num (22-27)
  savechar = buffer_[27];
  buffer_[27] = '\0';
  int resnum = atoi( buffer_+22 );
  buffer_[27] = savechar;
  return Residue(resnum, resname);
}

// PDBtype::pdb_XYZ()
// NOTE: Should check for NULL Xout
void PDBtype::pdb_XYZ(double *Xout) {
  // X coord (30-38)
  char savechar = buffer_[38];
  buffer_[38] = '\0';
  Xout[0] = atof( buffer_+30 );
  buffer_[38] = savechar;
  // Y coord (38-46)
  savechar = buffer_[46];
  buffer_[46] = '\0';
  Xout[1] = atof( buffer_+38 );
  buffer_[46] = savechar;
  // Z coord (46-54)
  savechar = buffer_[54];
  buffer_[54] = '\0';
  Xout[2] = atof( buffer_+46 );
  buffer_[54] = savechar;
}

// PDBtype::XYZ()
const double *PDBtype::XYZ() {
  // X coord (30-38)
  char savechar = buffer_[38];
  buffer_[38] = '\0';
  XYZ_[0] = atof( buffer_+30 );
  buffer_[38] = savechar;
  // Y coord (38-46)
  savechar = buffer_[46];
  buffer_[46] = '\0';
  XYZ_[1] = atof( buffer_+38 );
  buffer_[46] = savechar;
  // Z coord (46-54)
  savechar = buffer_[54];
  buffer_[54] = '\0';
  XYZ_[2] = atof( buffer_+46 );
  buffer_[54] = savechar;
  return XYZ_;
}

void PDBtype::pdb_write_ATOM(FileIO *IO, PDB_RECTYPE Record, int anum, NameType name,
                             NameType resnameIn, char chain, int resnum,
                             double X, double Y, double Z)
{
  pdb_write_ATOM(IO, Record, anum, name, resnameIn, chain, resnum, X, Y, Z, 0, 0, "", false);
}

/// Write out an ATOM or HETATM record
/** \return the number of characters written */
void PDBtype::pdb_write_ATOM(FileIO *IO, PDB_RECTYPE Record, int anum, NameType name,
                             NameType resnameIn, char chain, int resnum,
                             double X, double Y, double Z,
                             float Occ, float B, const char* End, bool highPrecision) 
{
  char resName[5], atomName[5];

  resName[4]='\0';
  atomName[4]='\0';
  // Residue number in PDB format can only be 4 digits wide
  while (resnum>9999) resnum-=9999;
  // Atom number in PDB format can only be 5 digits wide
  while (anum>99999) anum-=99999;
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
  // Atom names in PDB format start from col 14 when <= 3 chars, 13 when 4 chars.
  if (name[3]!=' ') { // 4 chars
    atomName[0] = name[0];
    atomName[1] = name[1];
    atomName[2] = name[2];
    atomName[3] = name[3];
  } else {            // <= 3 chars
    atomName[0] = ' ';
    atomName[1] = name[0];
    atomName[2] = name[1];
    atomName[3] = name[2];
  }

  char *ptr = buffer_;
  sprintf(ptr,"%-6s%5i %-4s%4s %c%4i",PDB_RECNAME[Record], anum, atomName, resName, chain, resnum);
  ptr+=26;
  if (Record == PDBTER) 
    sprintf(ptr,"\n");
  else if (highPrecision)
    sprintf(ptr,"    %8.3f%8.3f%8.3f%8.4f%8.4f%10s\n", X, Y, Z, Occ, B, End);
  else
    sprintf(ptr,"    %8.3f%8.3f%8.3f%6.2f%6.2f%14s\n", X, Y, Z, Occ, B, End);
  IO->Write(buffer_, strlen(buffer_));
}
 
/*/// PDB Record Chain ID
char PDB::pdb_chainID(char *buffer) {
  return buffer[21];
}*/

/*/// PDB Record Occupancy 
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
}*/

/// PDB Record B-factor
/*double pdb_Bfactor(char *buffer) {
  char temp[7];

  temp[0]=buffer[60];
  temp[1]=buffer[61];
  temp[2]=buffer[62];
  temp[3]=buffer[63];
  temp[4]=buffer[64];
  temp[5]=buffer[65];
  temp[6]='\0';
  return atof(temp);
}*/

/// 10 chars between Bfactor and element
/*char *pdb_lastChar(char *buffer) {
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
}*/

/// Element. If blank, try to guess from the name
/** Currently recognized: H C N O F Cl Br S P */
/*char *pdb_elt(char *buffer) {
  char *E,*ptr;
  char name[5];

  E=(char*) malloc(3*sizeof(char));
  E[0]=buffer[76]; // Element
  E[1]=buffer[77]; // Element
  E[2]='\0';*/

/*  if (E[0]==' ' && E[1]==' ') {
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
  }*/

/*  return E;
}*/

/// Charge. 
/*char *pdb_charge(char *buffer) {
  char *E;

  E=(char*) malloc(3*sizeof(char));
  E[0]=buffer[78]; // Charge
  E[1]=buffer[79]; // Charge
  E[2]='\0';
  return E;
}*/

