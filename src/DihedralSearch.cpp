#include "DihedralSearch.h"
#include "CpptrajStdio.h"

struct DIH_TYPE {
  int off;
  const char* an0;
  const char* an1;
  const char* an2;
  const char* an3;
  const char* name;
};

/** Recognized dihedral types go here. Must correspond to DihedralType,
  * except for NDIHTYPE which does not get an entry.
  */
static const DIH_TYPE DIH[] = {
  {-1, "C"  , "N"  , "CA" , "C"  , "phi"    }, // PHI: C0-N1-CA1-C1
  { 1, "N"  , "CA" , "C"  , "N"  , "psi"    }, // PSI: N0-CA0-C0-N1
  { 0, "N"  , "CA" , "CB" , "CG",  "chip"   }, // Protein CHI: N-CA-CB-CG
  {-2, "CA" , "C"  , "N"  , "CA" , "omega"  }, // OMEGA: CA0-C0-N1-CA1
  {-1, "O3'", "P"  , "O5'", "C5'", "alpha"  }, // ALPHA: 
  { 0, "P"  , "O5'", "C5'", "C4'", "beta"   }, // BETA:
  { 0, "O5'", "C5'", "C4'", "C3'", "gamma"  }, // GAMMA:
  { 0, "C5'", "C4'", "C3'", "O3'", "delta"  }, // DELTA:
  { 1, "C4'", "C3'", "O3'", "P"  , "epsilon"}, // EPSILON:
  { 2, "C3'", "O3'", "P"  , "O5'", "zeta"   }, // ZETA:
  { 0, "O4'", "C1'", "C2'", "C3'", "nu1"    }, // NU1:
  { 0, "C1'", "C2'", "C3'", "C4'", "nu2"    }, // NU2:
  { 0, "O4'", "C1'", "N9" , "C4",  "chin"   }  // Nucleic CHI:
};

void DihedralSearch::ListKnownTypes() {
  for (int i = 0; i < (int)NDIHTYPE; ++i)
    mprintf(" %s", DIH[i].name);
  mprintf("\n");
}

void DihedralSearch::OffsetHelp() {
  mprintf("\t\tOffset -2=<at0><at1> in previous res, -1=<at0> in previous res,\n");
  mprintf("\t\t        0=All <atX> in single res,\n");
  mprintf("\t\t        1=<at3> in next res, 2=<at2><at3> in next res.\n");
}

DihedralSearch::DihedralType DihedralSearch::GetType(std::string const& typeIn) {
  for (int i = 0; i < (int)NDIHTYPE; ++i)
    if (typeIn.compare(DIH[i].name)==0)
      return (DihedralType)i;
  return NDIHTYPE;
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask() : a0_(-1), a1_(-1), a2_(-1), 
                                               a3_(-1), res_(-1) {}

// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask(int a0, int a1, int a2, int a3, 
                                           int res, std::string const& name) :
  a0_(a0), a1_(a1), a2_(a2), a3_(a3), res_(res), name_(name) {}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - DihedralToken
DihedralSearch::DihedralToken::DihedralToken(int off, 
                                             NameType const& an0, NameType const& an1,
                                             NameType const& an2, NameType const& an3,
                                             std::string const& name) :
  offset_(off), aname0_(an0), aname1_(an1), aname2_(an2), aname3_(an3),
  name_(name)
{}

// DihedralSearch::DihedralToken::FindDihedralAtoms()
DihedralSearch::DihedralMask 
  DihedralSearch::DihedralToken::FindDihedralAtoms(Topology const& topIn, int resIn)
{
  int a1Res = resIn;
  int a2Res = resIn;
  int a3Res = resIn;
  int a4Res = resIn;
  switch (offset_) {
    case -2: --a2Res;        // -1 a2 and a1
    case -1: --a1Res; break; // -1 a1 only
    case  2: ++a3Res;        // -1 a3 and a4
    case  1: ++a4Res; break; // -1 a4 only
  }
  int atom1 = topIn.FindAtomInResidue(a1Res, aname0_);
  if (atom1 == -1) return DihedralMask();
  int atom2 = topIn.FindAtomInResidue(a2Res, aname1_);
  if (atom2 == -1) return DihedralMask();
  int atom3 = topIn.FindAtomInResidue(a3Res, aname2_);
  if (atom3 == -1) return DihedralMask();
  int atom4 = topIn.FindAtomInResidue(a4Res, aname3_);
  if (atom4 == -1) return DihedralMask();
  // All atoms found at this point.
  return DihedralMask(atom1, atom2, atom3, atom4, resIn, name_);
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
DihedralSearch::DihedralSearch() {}

// DihedralSearch::SearchFor()
int DihedralSearch::SearchFor(DihedralType typeIn) {
  dihedralTokens_.push_back( DihedralToken(DIH[typeIn].off, 
                               DIH[typeIn].an0, DIH[typeIn].an1, 
                               DIH[typeIn].an2, DIH[typeIn].an3,
                               DIH[typeIn].name) );
  return 0;
}

// DihedralSearch::SearchForArgs()
/** See if ArgList has any recognized dihedral type keywords. */
void DihedralSearch::SearchForArgs(ArgList& argIn) {
  for (int i = 0; i < (int)NDIHTYPE; ++i) {
    if (argIn.hasKey( DIH[i].name ))
      SearchFor( (DihedralType)i );
  }
}

// DihedralSearch::SearchForNewType()
/** Add new type to search for. */
int DihedralSearch::SearchForNewType(int off, std::string const& an0, std::string const& an1,
                                     std::string const& an2, std::string const& an3,
                                     std::string const& name)
{
  for (std::vector<DihedralToken>::iterator dih = dihedralTokens_.begin();
                                            dih != dihedralTokens_.end(); ++dih)
    if ( (*dih).Name() == name ) {
      mprintf("Warning: Dihedral type %s already defined.\n", name.c_str());
      return 1;
    }
  dihedralTokens_.push_back( DihedralToken(off, an0, an1, an2, an3, name) );
  return 0;
}

// DihedralSearch::SearchForAll()
/** If no dihedrals selected yet, select all. */
int DihedralSearch::SearchForAll() {
  if (!dihedralTokens_.empty()) return 0;
  for (int dih=0; dih < (int)NDIHTYPE; dih++)
    SearchFor((DihedralType)dih);
  return 0;
}

// DihedralSearch::FindDihedrals()
int DihedralSearch::FindDihedrals(Topology const& currentParm, Range const& rangeIn)
{
  dihedrals_.clear();
  for (Range::const_iterator res = rangeIn.begin(); res != rangeIn.end(); ++res)
  {
    for (std::vector<DihedralToken>::iterator dih = dihedralTokens_.begin();
                                              dih != dihedralTokens_.end(); ++dih)
    {
      dihedrals_.push_back( (*dih).FindDihedralAtoms(currentParm, *res) );
      if (dihedrals_.back().None()) {
        mprintf("Warning: Dihedral %s not found for residue %i\n", 
                (*dih).Name().c_str(), *res + 1);
        dihedrals_.pop_back();
      } 
    }
  }
  if (dihedrals_.empty()) {
    mprintf("Warning: No dihedrals selected for topology %s\n", currentParm.c_str());
    return 1;
  }
  //mprintf("\tFound %u dihedrals.\n", dihedrals_.size()); // DEBUG
  return 0;
}

// DihedralSearch::Clear()
void DihedralSearch::Clear() {
  dihedralTokens_.clear();
  dihedrals_.clear();
}

// DihedralSearch::ClearFound()
void DihedralSearch::ClearFound() {
  dihedrals_.clear();
}

// DihedralSearch::PrintTypes()
void DihedralSearch::PrintTypes() {
  for (std::vector<DihedralToken>::iterator dih = dihedralTokens_.begin();
                                            dih != dihedralTokens_.end(); ++dih)
  {
    mprintf(" %s", (*dih).Name().c_str());
  }
}

// VisitAtom()
static void VisitAtom( Topology const& topIn, int atm, std::vector<bool>& Visited )
{
  // If this atom has already been visited return
  if (Visited[atm]) return;
  // Mark this atom as visited
  Visited[atm] = true;
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = topIn[atm].bondbegin();
                           bondedatom != topIn[atm].bondend(); ++bondedatom)
    VisitAtom(topIn, *bondedatom, Visited);
}

// DihedralSearch::MovingAtoms()
AtomMask DihedralSearch::MovingAtoms(Topology const& topIn, int atom0, int atom1) {
  std::vector<bool> Visited( topIn.Natom(), false );
  // Mark atom0 as already visited
  Visited[atom0] = true;
  for (Atom::bond_iterator bndatm = topIn[atom1].bondbegin();
                           bndatm != topIn[atom1].bondend(); ++bndatm)
  {
    if ( *bndatm != atom0 )
      VisitAtom( topIn, *bndatm, Visited );
  }
  // Everything marked T will move.
  AtomMask Rmask;
  for (int maskatom = 0; maskatom < (int)Visited.size(); maskatom++) {
    if (Visited[maskatom])
      Rmask.AddAtom(maskatom);
  }
  return Rmask;
}
