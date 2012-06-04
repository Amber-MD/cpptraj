#include <cstring>
#include <cstdio>
#include "Mol2File.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Mol2File::Mol2File() : 
  mol2debug_(0),
  mol2atoms_(0),
  mol2bonds_(0)
{
  XYZ_[0] = 0;
  XYZ_[1] = 0;
  XYZ_[2] = 0;
}

/// Tripos Tags - must be in same order as enum type TRIPOSTAG
const char Mol2File::TRIPOSTAGTEXT[4][22]={
  "@<TRIPOS>MOLECULE",
  "@<TRIPOS>ATOM",
  "@<TRIPOS>BOND",
  "@<TRIPOS>SUBSTRUCTURE"
};

// Mol2File::IsMol2Keyword()
bool Mol2File::IsMol2Keyword() {
  if (strncmp(buffer_,"@<TRIPOS>",9)==0)
    return true;
  return false;
}

// Mol2File::GetLine()
bool Mol2File::GetLine(FileIO *IO) {
  return ( IO->Gets(buffer_, BUF_SIZE_) != 0 );
}

// Mol2File::ID()
bool Mol2File::ID( FileIO *IO ) {
  for (int line = 0; line < 10; line++) {
    if ( GetLine( IO ) ) return false;
    //mprintf("DEBUG: MOL2LINE %i: [%s]\n",line,buffer_);
    if ( IsMol2Keyword() ) {
      return true;
    }
  }
  return false;
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

// Mol2File::ReadMolecule()
/** Read in MOLECULE section of mol2file. Set title, number of atoms,
  * and number of bonds.
  */
bool Mol2File::ReadMolecule( FileIO *IO ) {
  // Scan to the section
  if ( ScanTo( IO, MOLECULE ) == 1 ) return true;
  //   Scan title
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return true;
  mol2title_.assign( buffer_ );
  RemoveTrailingWhitespace( mol2title_ );
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

// Mol2File::Mol2Bond()
void Mol2File::Mol2Bond(int &at1, int &at2) {
  sscanf(buffer_,"%*i %i %i\n", &at1, &at2);
}

// atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
Atom Mol2File::Mol2Atom() {
  char mol2name[10], mol2type[10];
  double mol2q;
  sscanf(buffer_, "%*i %s %*f %*f %*f %s %*i %*s %lf", mol2name, mol2type, &mol2q);
  NameType m2name( mol2name );
  // Replace all asterisks with single quote.
  m2name.ReplaceAsterisk();
  return Atom( m2name, mol2type, mol2q );
}

// Mol2File::Mol2Residue()
Residue Mol2File::Mol2Residue() {
  char resname[10];
  int resnum;
  sscanf(buffer_,"%*i %*s %*f %*f %*f %*s %i %s",&resnum,resname);
  NameType rname( resname );
  // Replace all asterisks with single quote.
  rname.ReplaceAsterisk();
  return Residue(resnum, rname);
}

// Mol2XYZ()
/** Given a Mol2 ATOM line, get the X Y and Z coords. */
void Mol2File::Mol2XYZ(double *X) {
  sscanf(buffer_,"%*i %*s %lf %lf %lf",X, X+1, X+2);
}

const double* Mol2File::XYZ() {
  sscanf(buffer_,"%*i %*s %lf %lf %lf",XYZ_, XYZ_+1, XYZ_+2);
  return XYZ_;
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

