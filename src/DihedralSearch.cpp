#include "DihedralSearch.h"
#include "CpptrajStdio.h"

const DihedralSearch::DihedralType DihedralSearch::D_FIRST = MetaData::ALPHA;
const DihedralSearch::DihedralType DihedralSearch::D_END   = MetaData::PUCKER;

/// Token to store pre-defined dihedral types.
struct DihedralSearch::DIH_TYPE {
  int off;
  DihedralType type;
  const char* an0;
  const char* an1;
  const char* an2;
  const char* an3;
};

/** Recognized dihedral types go here.  */
const DihedralSearch::DIH_TYPE DihedralSearch::DIH[] = {
  {-1, MetaData::PHI,     "C"  , "N"  , "CA" , "C"   }, // PHI: C0-N1-CA1-C1
  { 1, MetaData::PSI,     "N"  , "CA" , "C"  , "N"   }, // PSI: N0-CA0-C0-N1
  { 0, MetaData::CHIP,    "N"  , "CA" , "CB" , "CG"  }, // Protein CHI: R,N,D,Q,E,H,L,K,M,F,P,W,Y
  { 0, MetaData::CHIP,    "N"  , "CA" , "CB" , "SG"  }, // Protein CHI: C
  { 0, MetaData::CHIP,    "N"  , "CA" , "CB" , "CG1" }, // Protein CHI: I,V
  { 0, MetaData::CHIP,    "N"  , "CA" , "CB" , "OG"  }, // Protein CHI: S
  { 0, MetaData::CHIP,    "N"  , "CA" , "CB" , "OG1" }, // Protein CHI: T
  {-2, MetaData::OMEGA,   "CA" , "C"  , "N"  , "CA"  }, // OMEGA: CA0-C0-N1-CA1
  {-1, MetaData::ALPHA,   "O3'", "P"  , "O5'", "C5'" }, // ALPHA: 
  { 0, MetaData::BETA,    "P"  , "O5'", "C5'", "C4'" }, // BETA:
  { 0, MetaData::GAMMA,   "O5'", "C5'", "C4'", "C3'" }, // GAMMA:
  { 0, MetaData::DELTA,   "C5'", "C4'", "C3'", "O3'" }, // DELTA:
  { 1, MetaData::EPSILON, "C4'", "C3'", "O3'", "P"   }, // EPSILON:
  { 2, MetaData::ZETA,    "C3'", "O3'", "P"  , "O5'" }, // ZETA:
  { 0, MetaData::NU1,     "O4'", "C1'", "C2'", "C3'" }, // NU1: Nucleic pucker
  { 0, MetaData::NU2,     "C1'", "C2'", "C3'", "C4'" }, // NU2: Nucleic pucker
  { 0, MetaData::CHIN,    "O4'", "C1'", "N9",  "C4"  }, // Nucleic CHI: Purine (A, G)
  { 0, MetaData::CHIN,    "O4'", "C1'", "N1",  "C2"  }, // Nucleic CHI: Pyrimidine (U, T, C)
  { 0, MetaData::H1P,     "H1'", "C1'", "N9",  "C4"  }, // Nucleic H1' sugar pucker-base (purine)
  { 0, MetaData::H1P,     "H1'", "C1'", "N1",  "C2"  }, // Nucleic H1' sugar pucker-base (pyrim.)
  { 0, MetaData::C2P,     "C2'", "C1'", "N9",  "C4"  }, // Nucleic C2' sugar pucker-base (purine)
  { 0, MetaData::C2P,     "C2'", "C1'", "N1",  "C2"  }, // Nucleic C2' sugar pucker-base (pyrim.)
  { 0, MetaData::UNDEFINED,   0,     0,     0,     0 }
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
  offset_(off),
  name_(name),
  type_(MetaData::UNDEFINED)
{
  aname_[0] = an0;
  aname_[1] = an1;
  aname_[2] = an2;
  aname_[3] = an3;
}

// CONSTRUCTOR - Recognized type 
DihedralSearch::DihedralToken::DihedralToken(DIH_TYPE const& dih) :
  offset_(dih.off),
  name_(MetaData::TypeString(dih.type)),
  type_(dih.type)
{
  aname_[0] = dih.an0;
  aname_[1] = dih.an1;
  aname_[2] = dih.an2;
  aname_[3] = dih.an3;
}

static inline int FindNameBondedTo(Atom const& atm, Topology const& top,
                                    NameType const& aname, int rnum)
{
  int atomX = -1;
  for (Atom::bond_iterator bat = atm.bondbegin(); bat != atm.bondend(); ++bat)
    if (top[*bat].Name() == aname && top[*bat].ResNum() == rnum)
    {
      atomX = *bat;
      break;
    }
  return atomX;
}

// DihedralSearch::DihedralToken::FindDihedralAtoms()
// FIXME handle cyclic peptides etc
DihedralSearch::DihedralMask 
  DihedralSearch::DihedralToken::FindDihedralAtoms(Topology const& topIn, int resIn) const
{
  int resnum[4] = {resIn, resIn, resIn, resIn};
  switch (offset_) {
    case -2: --resnum[1];        // -1 a2 and a1
    case -1: --resnum[0]; break; // -1 a1 only
    case  2: ++resnum[2];        // +1 a3 and a4
    case  1: ++resnum[3]; break; // +1 a4 only
  }
  if (resnum[0] < 0 || resnum[0] >= topIn.Nres()) return DihedralMask();
  // Find the first atom in specified residue.
  int atom0 = -1;
  for (int at = topIn.Res(resnum[0]).FirstAtom(); at != topIn.Res(resnum[0]).LastAtom(); at++)
    if (topIn[at].Name() == aname_[0]) {
      atom0 = at;
      break;
    }
  if (atom0 < 0) return DihedralMask();
  // Subsequent atoms must be bonded to the previous atom.
  int atom1 = FindNameBondedTo(topIn[atom0], topIn, aname_[1], resnum[1]);
  if (atom1 < 0) return DihedralMask();
  int atom2 = FindNameBondedTo(topIn[atom1], topIn, aname_[2], resnum[2]);
  if (atom2 < 0) return DihedralMask();
  int atom3 = FindNameBondedTo(topIn[atom2], topIn, aname_[3], resnum[3]);
  if (atom3 < 0) return DihedralMask();
  // All atoms found at this point.
  return DihedralMask(atom0, atom1, atom2, atom3, resIn, name_, type_);
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
DihedralSearch::DihedralSearch() {}

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
