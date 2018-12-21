#include "DihedralSearch.h"
#include "CpptrajStdio.h"

const DihedralSearch::DihedralType DihedralSearch::D_FIRST = MetaData::ALPHA;
const DihedralSearch::DihedralType DihedralSearch::D_END   = MetaData::PUCKER;

/// Token to store pre-defined dihedral types.
struct DihedralSearch::DIH_TYPE {
  int cidx;          ///< Index of the "central" atom (i.e. what res should it belong to)
  DihedralType type; ///< Dihedral MetaData type
  const char* an0;   ///< First atom name
  const char* an1;   ///< Second atom name
  const char* an2;   ///< Third atom name
  const char* an3;   ///< Fourth atom name
};

/** Recognized dihedral types go here.  */
const DihedralSearch::DIH_TYPE DihedralSearch::DIH[] = {
  { 2, MetaData::PHI,     "C"  , "N"  , "CA" , "C"   }, // PHI: C0-N1-CA1-C1
  { 2, MetaData::PSI,     "N"  , "CA" , "C"  , "N"   }, // PSI: N0-CA0-C0-N1
  { 2, MetaData::CHIP,    "N"  , "CA" , "CB" , "CG"  }, // Protein CHI: R,N,D,Q,E,H,L,K,M,F,P,W,Y
  { 2, MetaData::CHIP,    "N"  , "CA" , "CB" , "SG"  }, // Protein CHI: C
  { 2, MetaData::CHIP,    "N"  , "CA" , "CB" , "CG1" }, // Protein CHI: I,V
  { 2, MetaData::CHIP,    "N"  , "CA" , "CB" , "OG"  }, // Protein CHI: S
  { 2, MetaData::CHIP,    "N"  , "CA" , "CB" , "OG1" }, // Protein CHI: T
  { 2, MetaData::CHI2,    "CA" , "CB" , "CG" , "CD"  }, // Protein CHI2: R, Q, E, I, K, P
  { 2, MetaData::CHI2,    "CA" , "CB" , "CG" , "OD1" }, // Protein CHI2: N, D
  { 2, MetaData::CHI2,    "CA" , "CB" , "CG" , "ND1" }, // Protein CHI2: H
  { 2, MetaData::CHI2,    "CA" , "CB" , "CG" , "CD1" }, // Protein CHI2: L, F, W, Y
  { 2, MetaData::CHI2,    "CA" , "CB" , "CG" , "SD"  }, // Protein CHI2: M
  { 2, MetaData::CHI3,    "CB" , "CG" , "CD" , "NE"  }, // Protein CHI3: R
  { 2, MetaData::CHI3,    "CB" , "CG" , "CD" , "OE1" }, // Protein CHI3: Q, E
  { 2, MetaData::CHI3,    "CB" , "CG" , "CD" , "CE"  }, // Protein CHI3: K
  { 2, MetaData::CHI3,    "CB" , "CG" , "SD" , "CE"  }, // Protein CHI3: M
  { 2, MetaData::CHI4,    "CG" , "CD" , "NE" , "CZ"  }, // Protein CHI4: R
  { 2, MetaData::CHI4,    "CG" , "CD" , "CE" , "NZ"  }, // Protein CHI4: K
  { 2, MetaData::CHI5,    "CD" , "NE" , "CZ" , "NH1" }, // Protein CHI5: R
  { 2, MetaData::OMEGA,   "CA" , "C"  , "N"  , "CA"  }, // OMEGA: CA0-C0-N1-CA1
  { 2, MetaData::ALPHA,   "O3'", "P"  , "O5'", "C5'" }, // ALPHA:
  { 2, MetaData::BETA,    "P"  , "O5'", "C5'", "C4'" }, // BETA:
  { 2, MetaData::GAMMA,   "O5'", "C5'", "C4'", "C3'" }, // GAMMA:
  { 2, MetaData::DELTA,   "C5'", "C4'", "C3'", "O3'" }, // DELTA:
  { 2, MetaData::EPSILON, "C4'", "C3'", "O3'", "P"   }, // EPSILON:
  { 1, MetaData::ZETA,    "C3'", "O3'", "P"  , "O5'" }, // ZETA:
  { 2, MetaData::NU0,     "C4'", "O4'", "C1'", "C2'" }, // NU0: Nucleic pucker
  { 2, MetaData::NU1,     "O4'", "C1'", "C2'", "C3'" }, // NU1: Nucleic pucker
  { 2, MetaData::NU2,     "C1'", "C2'", "C3'", "C4'" }, // NU2: Nucleic pucker
  { 2, MetaData::NU3,     "C2'", "C3'", "C4'", "O4'" }, // NU3: Nucleic pucker
  { 2, MetaData::NU4,     "C3'", "C4'", "O4'", "C1'" }, // NU4: Nucleic pucker
  { 2, MetaData::CHIN,    "O4'", "C1'", "N9",  "C4"  }, // Nucleic CHI: Purine (A, G)
  { 2, MetaData::CHIN,    "O4'", "C1'", "N1",  "C2"  }, // Nucleic CHI: Pyrimidine (U, T, C)
  { 2, MetaData::H1P,     "H1'", "C1'", "N9",  "C4"  }, // Nucleic H1' sugar pucker-base (purine)
  { 2, MetaData::H1P,     "H1'", "C1'", "N1",  "C2"  }, // Nucleic H1' sugar pucker-base (pyrim.)
  { 2, MetaData::C2P,     "C2'", "C1'", "N9",  "C4"  }, // Nucleic C2' sugar pucker-base (purine)
  { 2, MetaData::C2P,     "C2'", "C1'", "N1",  "C2"  }, // Nucleic C2' sugar pucker-base (pyrim.)
  { 2, MetaData::UNDEFINED,   0,     0,     0,     0 }
};

// DihedralSearch::ListKnownTypes()
void DihedralSearch::ListKnownTypes() {
  for (int i = (int)D_FIRST; i < (int)D_END; ++i)
    mprintf(" %s", MetaData::TypeString((DihedralType)i));
  mprintf("\n");
}

// TODO: const char*
void DihedralSearch::OffsetHelp() {
  mprintf("\t\tOffset -2=<a0><a1> in previous res, -1=<a0> in previous res,\n"
          "\t\t        0=All <aX> in single res,\n"
          "\t\t        1=<a3> in next res, 2=<a2><a3> in next res.\n");
}

// DihedralSearch::GetType()
DihedralSearch::DihedralType DihedralSearch::GetType(std::string const& typeIn) {
  for (int i = (int)D_FIRST; i < (int)D_END; ++i)
    if (typeIn.compare( MetaData::TypeString((DihedralType)i) )==0)
      return (DihedralType)i;
  return MetaData::UNDEFINED;
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask() : 
  a0_(-1), a1_(-1), a2_(-1), a3_(-1), res_(-1), type_(MetaData::UNDEFINED) {}

// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask(int a0, int a1, int a2, int a3, 
                                           int res, std::string const& n,
                                           DihedralType t) :
  a0_(a0), a1_(a1), a2_(a2), a3_(a3), res_(res), name_(n), type_(t) {}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - Custom type 
DihedralSearch::DihedralToken::DihedralToken(int off, 
                                             NameType const& an0, NameType const& an1,
                                             NameType const& an2, NameType const& an3,
                                             std::string const& name) :
  centerIdx_(2),
  name_(name),
  type_(MetaData::UNDEFINED)
{
  // If offset is 2, central atom is second, otherwise use third.
  if (off > 1) centerIdx_ = 1;
  aname_[0] = an0;
  aname_[1] = an1;
  aname_[2] = an2;
  aname_[3] = an3;
}

// CONSTRUCTOR - Recognized type 
DihedralSearch::DihedralToken::DihedralToken(DIH_TYPE const& dih) :
  centerIdx_(dih.cidx),
  name_(MetaData::TypeString(dih.type)),
  type_(dih.type)
{
  aname_[0] = dih.an0;
  aname_[1] = dih.an1;
  aname_[2] = dih.an2;
  aname_[3] = dih.an3;
}

/// \return index of atom named aname bonded to atm in Topology top, or -1 if not found.
static inline int FindNameBondedTo(Atom const& atm, Topology const& top, NameType const& aname)
{
  int atomX = -1;
  for (Atom::bond_iterator bat = atm.bondbegin(); bat != atm.bondend(); ++bat)
    if (top[*bat].Name() == aname)
    {
      atomX = *bat;
      break;
    }
  return atomX;
}

// DihedralSearch::DihedralToken::FindDihedralAtoms()
/** Search given residue for this dihedral.
  * \return Appropriate DihedralMask if found, or empty DihedralMask if not found.
  */
DihedralSearch::DihedralMask 
  DihedralSearch::DihedralToken::FindDihedralAtoms(Topology const& topIn, int resIn) const
{
  Residue const& Res = topIn.Res(resIn);
  int atIdx[4];
  // See if central atom exists.
  atIdx[centerIdx_] = -1;
  for (int at = Res.FirstAtom(); at != Res.LastAtom(); at++)
  {
    if (topIn[at].Name() == aname_[centerIdx_]) {
      atIdx[centerIdx_] = at;
      break;
    }
  }
  if (atIdx[centerIdx_] == -1) return DihedralMask();
  // Find atoms in forward direction.
  for (int i = centerIdx_+1; i < 4; i++)
  {
    int idx = FindNameBondedTo(topIn[atIdx[i-1]], topIn, aname_[i]);
    if (idx == -1) return DihedralMask();
    atIdx[i] = idx;
  }
  // Find atoms in reverse direction.
  for (int i = centerIdx_-1; i > -1; i--)
  {
    int idx = FindNameBondedTo(topIn[atIdx[i+1]], topIn, aname_[i]);
    if (idx == -1) return DihedralMask();
    atIdx[i] = idx;
  }
  // All atoms found at this point.
  return DihedralMask(atIdx[0], atIdx[1], atIdx[2], atIdx[3], resIn, name_, type_);
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
DihedralSearch::DihedralSearch() {}

/** COPY CONSTRUCTOR - Set up for same tokens. */
DihedralSearch::DihedralSearch(DihedralSearch const& rhs) :
  dihedralTokens_(rhs.dihedralTokens_)
{}

// DihedralSearch::SearchFor()
/** Search for all types matching typeIn. */
int DihedralSearch::SearchFor(DihedralType typeIn) {
  for (DIH_TYPE const* ptr = DIH; ptr->type != MetaData::UNDEFINED; ++ptr)
    if (ptr->type == typeIn)
      dihedralTokens_.push_back( DihedralToken( *ptr ));
  return 0;
}

// DihedralSearch::SearchForArgs()
/** See if ArgList has any recognized dihedral type keywords. */
void DihedralSearch::SearchForArgs(ArgList& argIn) {
  for (int i = (int)D_FIRST; i < (int)D_END; ++i) {
    if (argIn.hasKey( MetaData::TypeString((DihedralType)i) ))
      SearchFor( (DihedralType)i );
  }
}

// DihedralSearch::SearchForNewType()
/** Add new type to search for. */
int DihedralSearch::SearchForNewType(int off, std::string const& an0, std::string const& an1,
                                     std::string const& an2, std::string const& an3,
                                     std::string const& name)
{
  for (std::vector<DihedralToken>::iterator tkn = dihedralTokens_.begin();
                                            tkn != dihedralTokens_.end(); ++tkn)
    if ( tkn->Name() == name ) {
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
  for (int dih = (int)D_FIRST; dih < (int)D_END; ++dih)
    SearchFor((DihedralType)dih);
  return 0;
}

// DihedralSearch::FindDihedrals()
int DihedralSearch::FindDihedrals(Topology const& currentParm, Range const& rangeIn)
{
  dihedrals_.clear();
  for (Range::const_iterator res = rangeIn.begin(); res != rangeIn.end(); ++res)
  {
    //std::string notFound; //TODO record this 
    for (std::vector<DihedralToken>::const_iterator tkn = dihedralTokens_.begin();
                                                    tkn != dihedralTokens_.end(); ++tkn)
    {
      dihedrals_.push_back( tkn->FindDihedralAtoms(currentParm, *res) );
      if (dihedrals_.back().None()) {
        //notFound.append( " " + tkn->Name() );
        dihedrals_.pop_back();
      } 
    }
    //if (!notFound.empty())
    //  mprintf("Warning: Dihedral%s not found for residue %i\n", notFound.c_str(), *res + 1);
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

// DihedralSearch::PrintTypes()
void DihedralSearch::PrintTypes() {
  for (std::vector<DihedralToken>::iterator tkn = dihedralTokens_.begin();
                                            tkn != dihedralTokens_.end(); ++tkn)
    mprintf(" %s", tkn->Name().c_str());
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
  std::vector<int> Rmask;
  for (int maskatom = 0; maskatom < (int)Visited.size(); maskatom++) {
    if (Visited[maskatom])
      Rmask.push_back(maskatom);
  }
  return AtomMask(Rmask, topIn.Natom());
}
