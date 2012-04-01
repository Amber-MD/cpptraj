#include <cstring>
#include <cstdio>
#include "Mol2File.h"
#include "CpptrajStdio.h"

Mol2File::Mol2File() : 
  mol2debug_(0),
  mol2atoms_(0),
  mol2bonds_(0)
{}

/// Tripos Tags - must be in same order as enum type TRIPOSTAG
const char Mol2File::TRIPOSTAGTEXT[4][22]={
  "@<TRIPOS>MOLECULE",
  "@<TRIPOS>ATOM",
  "@<TRIPOS>BOND",
  "@<TRIPOS>SUBSTRUCTURE"
};

bool Mol2File::IsMol2Keyword() {
  if (strncmp(buffer_,"@<TRIPOS>",9)==0)
    return true;
  return false;
}

bool Mol2File::GetLine(FileIO *IO) {
  return ( IO->Gets(buffer_, BUF_SIZE_) == 0 );
}

// Mol2File::ScanTo()
/** Scan to the specified TRIPOS section of file.
  * \return 0 if the tag was found, 1 if not found.
  */
int Mol2File::ScanTo( FileIO *IO, TRIPOSTAG tag ) {
  int tagSize = (int)strlen(TRIPOSTAGTEXT[tag]);
  while ( IO->Gets(buffer_,BUF_SIZE_)==0 ) {
    //mprintf("DEBUG: Line [%s]\n",buffer);
    //mprintf("DEBUG: Targ [%s]\n",TRIPOSTAGTEXT[tag]); 
    if (strncmp(buffer_,TRIPOSTAGTEXT[tag],tagSize)==0) return 0;
  }
  // Suppress this warning so routine can be used to scan # frames
  //mprintf("Warning: Mol2File::ScanTo(): Could not find tag %s\n",TRIPOSTAGTEXT[tag]);
  return 1;
}

bool Mol2File::ReadMolecule( FileIO *IO ) {
  // Scan to the section
  if ( ScanTo( IO, MOLECULE ) == 1 ) return true;
  //   Scan title
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return true;
  mol2title_.assign( buffer_ );
  if (mol2debug_>0) mprintf("      Mol2 Title: [%s]\n",mol2title_.c_str());
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return true;
  mol2atoms_ = 0;
  mol2bonds_ = 0;
  sscanf(buffer_,"%i %i",&mol2atoms_, &mol2bonds_);
  if (mol2debug_>0) {
    mprintf("      Mol2 #atoms: %i\n",mol2atoms_);
    mprintf("      Mol2 #bonds: %i\n",mol2bonds_);
  }
  return false;
}

/** Used to only read # atoms in next MOLECULE record. 
  * \return # atoms in next MOLECULE, -1 on error or end of file.
  */
int Mol2File::NextMolecule( FileIO *IO ) {
  int natom;
  // Scan to the section
  if ( ScanTo( IO, MOLECULE ) == 1 ) return -1;
  // Scan past the title
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return -1;
  // Scan # atoms
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return -1;
  sscanf(buffer_, "%i", &natom);
  return natom;
}

void Mol2File::Mol2Bond(int &at1, int &at2) {
  sscanf(buffer_,"%*i %i %i\n", &at1, &at2);
}

// atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
Atom Mol2File::Mol2Atom() {
  char mol2name[10], mol2type[10];
  double mol2q;
  double XYZ_[3];
  sscanf(buffer_, "%*i %s %lf %lf %lf %s %*i %*s %lf", mol2name,XYZ_, XYZ_+1, XYZ_+2, 
         mol2type, &mol2q);
  return Atom( mol2name, XYZ_, mol2type, mol2q );
}

Residue Mol2File::Mol2Residue() {
  char resname[10];
  int resnum;
  sscanf(buffer_,"%*i %*s %*f %*f %*f %*s %i %s",&resnum,resname);
  return Residue(resnum, resname);
}

// Mol2XYZ()
/** Given a Mol2 ATOM line, get the X Y and Z coords. */
void Mol2File::Mol2XYZ(double *X) {
  sscanf(buffer_,"%*i %*s %lf %lf %lf",X, X+1, X+2);
}

void Mol2File::SetMol2Natoms(int natomIn) {
  mol2atoms_ = natomIn;
}

void Mol2File::SetMol2Nbonds(int nbondIn) {
  mol2bonds_ = nbondIn;
}

int Mol2File::Mol2Natoms() {
  return mol2atoms_;
}

int Mol2File::Mol2Nbonds() {
  return mol2bonds_;
}

std::string& Mol2File::Mol2Title() {
  return mol2title_;
}

// Mol2AtomName()
/** Given a Mol2 ATOM line, return atom name. Trim to 4 chars to be consistent 
  * with the rest of Amber. 
  */
/*
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
}*/

// Mol2AtomType
/** Given a Mol2 ATOM line, return atom type. Try to convert Sybyl atom type
  * to amber type.
  * Sybyl atom types seem to have at most 5 chars - increasing size of NAME
  * huld be ok.
  */
/*
int Mol2AtomType(char *buffer, NAME type) {
  char tmp[10];
  if (buffer==NULL || type==NULL) return 1;
  sscanf(buffer,"%*i %*s %*f %*f %*f %s",tmp);
*/
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
/*
    // Copy over the first 5 chars
    strncpy(type,tmp,5);
    type[5]='\0';
  //}
  return 0;
}
*/
// Mol2ResNumName()
/** Get residue number and name, trim to 4 chars
  */
/*int Mol2ResNumName(char *buffer, int *resnum, NAME resname) {
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
}*/

// Mol2Charge()
/** Get charge */
/*double Mol2Charge(char *buffer) {
  double q;
  sscanf(buffer,"%*i %*s %*f %*f %*f %*s %*i %*s %lf",&q);
  return q;
} */

