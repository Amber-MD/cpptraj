#include <cstring> //strlen
#include <cstdio> //sscanf
#include "Mol2File.h"
#include "CpptrajStdio.h" // To print debug info
#include "StringRoutines.h" // RemoveTrailingWhitespace

// CONSTRUCTOR
Mol2File::Mol2File() : 
  mol2debug_(0),
  mol2atoms_(0),
  mol2bonds_(0)
{}

/// Tripos Tags - must be in same order as enum type TRIPOSTAG
const char* Mol2File::TRIPOSTAGTEXT[]={
  "@<TRIPOS>MOLECULE",
  "@<TRIPOS>ATOM",
  "@<TRIPOS>BOND",
  "@<TRIPOS>SUBSTRUCTURE"
};

// Mol2File::IsMol2Keyword()
bool Mol2File::IsMol2Keyword(const char* key) {
  if (strncmp(key, "@<TRIPOS>",9)==0)
    return true;
  return false;
}

// Mol2File::ID_Mol2()
bool Mol2File::ID_Mol2(CpptrajFile& fileIn) {
  // NOTE: ASSUMES FILE IS ALREADY SETUP!
  if (fileIn.OpenFile()) return false;
  for (int line = 0; line < 10; line++) {
    std::string nextLine = fileIn.GetLine();
    //mprintf("DEBUG: MOL2LINE %i: [%s]\n",line,linebuffer_);
    if ( IsMol2Keyword(nextLine.c_str()) ) {
      fileIn.CloseFile();
      return true;
    }
  }
  fileIn.CloseFile();
  return false;
}

// Mol2File::ScanTo()
/** \return 0 if the tag was found, 1 if not found. */
int Mol2File::ScanTo( TRIPOSTAG tag ) {
  int tagSize = (int)strlen(TRIPOSTAGTEXT[tag]);
  //mprintf("DEBUG: SCANNING TO MOL2 TAG '%s'\n", TRIPOSTAGTEXT[tag]);
  while ( Gets(linebuffer_, BUF_SIZE)==0 ) {
    //mprintf("DEBUG: Line [%s]\n",linebuffer_);
    if (strncmp(linebuffer_, TRIPOSTAGTEXT[tag], tagSize)==0) return 0;
  }
  // Suppress this warning so routine can be used to scan # frames
  //mprintf("Warning: Mol2File::ScanTo(): Could not find tag %s\n",TRIPOSTAGTEXT[tag]);
  return 1;
}

void Mol2File::WriteHeader( Mol2File::TRIPOSTAG tag ) {
  Printf("%s\n", TRIPOSTAGTEXT[tag] );
}

// Mol2File::ReadMolecule()
/** Set title, number of atoms, and number of bonds. */
bool Mol2File::ReadMolecule( ) {
  // Scan to the section
  if ( ScanTo( MOLECULE ) == 1 ) return true;
  //   Scan title
  if ( Gets(linebuffer_, BUF_SIZE) ) return true;
  mol2title_.assign( linebuffer_ );
  RemoveTrailingWhitespace( mol2title_ );
  if (mol2debug_>0) mprintf("      Mol2 Title: [%s]\n",mol2title_.c_str());
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( Gets(linebuffer_, BUF_SIZE) ) return true;
  mol2atoms_ = 0;
  mol2bonds_ = 0;
  if (sscanf(linebuffer_,"%i %i",&mol2atoms_, &mol2bonds_) != 2) {
    mprinterr("Error: Mol2File: Could not read # atoms/ # bonds.\n");
    return false;
  }
  if (mol2debug_>0) {
    mprintf("\tMol2 #atoms: %i\n",mol2atoms_);
    mprintf("\tMol2 #bonds: %i\n",mol2bonds_);
  }
  return false;
}

bool Mol2File::WriteMolecule(bool hasCharges, int mol2res) {
  Printf("%s\n", TRIPOSTAGTEXT[MOLECULE]);
  // mol_name
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  // mol_type
  // charge_type
  // [status_bits
  // [mol_comment]]
  Printf("%s\n", mol2title_.c_str());
  Printf("%5i %5i %5i %5i %5i\n", mol2atoms_, mol2bonds_, mol2res, 0, 0);
  Printf("SMALL\n"); // May change this later
  if ( hasCharges )
    Printf("USER_CHARGES\n"); // May change this later
  else
    Printf("NO_CHARGES\n");
  Printf("\n\n");
  return false;
}

/** \return # atoms in next MOLECULE, -1 on error or end of file. */
int Mol2File::NextMolecule( ) {
  int natom = 0;
  // Scan to the section
  if ( ScanTo( MOLECULE ) == 1 ) return -1;
  // Scan past the title
  if ( Gets(linebuffer_, BUF_SIZE) ) return -1;
  // Scan # atoms
  if ( Gets(linebuffer_, BUF_SIZE) ) return -1;
  sscanf(linebuffer_, "%i", &natom);
  return natom;
}

// Mol2File::Mol2Bond()
int Mol2File::Mol2Bond(int& at1, int& at2) {
  if ( Gets(linebuffer_, BUF_SIZE) != 0 ) return 1;
  // bond_id origin_atom_id target_atom_id bond_type [status_bits]
  sscanf(linebuffer_,"%*i %i %i\n", &at1, &at2);
  return 0;
}

// Mol2File::Mol2XYZ()
int Mol2File::Mol2XYZ(double *X) {
  if ( Gets(linebuffer_, BUF_SIZE) != 0 ) return 1;
  sscanf(linebuffer_,"%*i %*s %lf %lf %lf",X, X+1, X+2);
  return 0;
}

// Mol2File::Mol2Atom()
Atom Mol2File::Mol2Atom() {
  char mol2name[10], mol2type[10];
  double mol2q;
  // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
  sscanf(linebuffer_, "%*i %s %*f %*f %*f %s %*i %*s %lf", mol2name, mol2type, &mol2q);
  NameType m2name( mol2name );
  // Replace all asterisks with single quote.
  m2name.ReplaceAsterisk();
  return Atom( m2name, mol2type, mol2q );
}

// Mol2File::Mol2Residue()
Residue Mol2File::Mol2Residue() {
  char resname[10];
  int current_res;
  sscanf(linebuffer_,"%*i %*s %*f %*f %*f %*s %i %s", &current_res, resname);
  NameType rname( resname );
  // Replace all asterisks with single quote.
  rname.ReplaceAsterisk();
  return Residue(rname, current_res, ' ', ' ');
}

void Mol2File::WriteMol2Atom(int atnum, Atom const& atomIn,
                             int resnum, const char* rname,
                             const double* Xptr)
{
  NameType atype;
  // If mapping is defined, try to use it.
  if (!Atype_to_Sybyl_.empty()) {
    AtypeMap::const_iterator it = Atype_to_Sybyl_.find( atomIn.Type() );
    if (it == Atype_to_Sybyl_.end()) {
      mprintf("Warning: SYBYL type for atom %i '%s' not found.\n", atnum, *atomIn.Type());
      atype = atomIn.Name();
    } else
      atype = it->second;
  } else {
    // If atom type is blank, set to atom name
    atype = atomIn.Type();
    if (atype == "")
      atype = atomIn.Name();
  }
  Printf("%7i %-8s %9.4lf %9.4lf %9.4lf %-5s %6i %-6s %10.6lf\n",
         atnum, atomIn.c_str(), Xptr[0], Xptr[1], Xptr[2],
         *atype, resnum, rname, atomIn.Charge());
}

const char* Mol2File::SYBYL_BOND_[] = {"1", "2", "3", "am", "ar"};

void Mol2File::WriteMol2Bond(int bnum, int at1, int at2,
                             NameType const& type1, NameType const& type2)
{
  int bidx = 0;
  // If mapping is defined, try to use it
  if (!Apair_to_Bond_.empty()) {
    // Always order the pair so first < second
    AtomPair bp;
    if (type1 < type2)
      bp = AtomPair(type1, type2);
    else
      bp = AtomPair(type2, type1);
    BndMap::const_iterator it = Apair_to_Bond_.find( bp );
    if (it != Apair_to_Bond_.end())
    //  mprintf("Warning: SYBYL bond for atom %i to %i (%s -- %s) not found.\n",
    //          at1+1, at2+1, *type1, *type2);
    //else
      bidx = it->second;
  }
  Printf("%5d %5d %5d %s\n", bnum, at1, at2, SYBYL_BOND_[bidx]);
}

void Mol2File::WriteMol2Substructure(int rnum, const char* rname, int firstatom) {
  Printf("%7d %4s %14d ****               0 ****  **** \n", rnum, rname, firstatom);
}

void Mol2File::ClearAmberMapping() {
  Atype_to_Sybyl_.clear();
  Apair_to_Bond_.clear();
}

int Mol2File::ReadAmberMapping(FileName const& AtomFile, FileName const& BondFile, int debug)
{
  CpptrajFile infile;
  const char* ptr = 0;
  char Atype[6]; // Amber Atom type: Max 2 letters
  char Stype[6]; // SYBYL Atom type: Max 5 letters
  if (!AtomFile.empty()) {
    // Expected format: <Amber Atom Type> <SYBYL Atom Type>
    if (infile.OpenRead( AtomFile )) return 1;
    ptr = infile.NextLine();
    std::pair<AtypeMap::iterator, bool> ret;
    typedef std::pair<NameType, NameType> Mpair;
    while (ptr != 0) {
      sscanf(ptr, "%s %s", Atype, Stype);
      NameType amber_type(Atype);
      NameType sybyl_type(Stype);
      ret = Atype_to_Sybyl_.insert( Mpair(amber_type, sybyl_type) );
      if (!ret.second) {
        if (sybyl_type != ret.first->second) {
          mprinterr("Error: Duplicate Amber atom type '%s' in '%s' has different SYBYL\n"
                    "Error:   has different SYBYL type '%s' than previous '%s'\n",
                    *amber_type, AtomFile.full(), *sybyl_type, *(ret.first->second));
          return 1;
        }
        mprintf("Warning: Duplicate Amber atom type '%s' in '%s'\n", *amber_type, AtomFile.full());
      }
      ptr = infile.NextLine();
    }
    infile.CloseFile();
    if (debug > 0) {
      mprintf("DEBUG: Atype_to_Sybyl has %zu values:\n", Atype_to_Sybyl_.size());
      for (AtypeMap::const_iterator ix = Atype_to_Sybyl_.begin(); ix != Atype_to_Sybyl_.end(); ++ix)
        mprintf("\t'%s' => '%s'\n", *(ix->first), *(ix->second));
    }
  }
  if (!BondFile.empty()) {
    // Expected format: <Amber Atom1 Type> <Amber Atom2 Type> <SYBYL Bond Type>
    if (infile.OpenRead( BondFile )) return 1;
    ptr = infile.NextLine();
    char Atype2[6]; // Amber Atom type: Max 2 letters
    std::pair<BndMap::iterator, bool> bret;
    typedef std::pair<AtomPair, int> Bpair;
    while (ptr != 0) {
      sscanf(ptr, "%s %s %s", Atype, Atype2, Stype);
      NameType a1(Atype);
      NameType a2(Atype2);
      // Always order the pair so first < second
      AtomPair bp;
      if (a1 < a2)
        bp = AtomPair(a1, a2);
      else
        bp = AtomPair(a2, a1);
      // Determine SYBYL bond type
      int btype = -1;
      if      (Stype[0] == '1') btype = 0;
      else if (Stype[0] == '2') btype = 1;
      else if (Stype[0] == '3') btype = 2;
      else if (Stype[0] == 'a') {
        if      (Stype[1] == 'm') btype = 3;
        else if (Stype[1] == 'r') btype = 4;
      }
      if (btype < 0) {
        mprinterr("Error: File '%s' contains unsupported SYBYL bond type '%s'\n",
                  BondFile.full(), Stype);
        return 1;
      }
      bret = Apair_to_Bond_.insert( Bpair(bp, btype) );
      if (!bret.second) {
        // Only make this an error if the type does not match.
        if (btype != bret.first->second) {
          mprinterr("Error: Duplicate bond '%s'-'%s' in '%s'\n"
                    "Error:   has different type %i than previous %i\n",
                    *a1, *a2, BondFile.full(), btype, bret.first->second);
          return 1;
        }
        mprintf("Warning: Duplicate bond '%s'-'%s' in '%s'\n", *a1, *a2, BondFile.full());
      }
      ptr = infile.NextLine();
    }
    infile.CloseFile();
    if (debug > 0) {
      mprintf("DEBUG: Apair_to_Bond has %zu values:\n", Apair_to_Bond_.size());
      for (BndMap::const_iterator ib = Apair_to_Bond_.begin(); ib != Apair_to_Bond_.end(); ++ib)
        mprintf("'%s'--'%s' => %s\n", *(ib->first.first), *(ib->first.second),
                SYBYL_BOND_[ib->second]);
    }
  }
  return 0;
}
