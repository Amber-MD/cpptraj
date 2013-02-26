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
  {-1, "C", "N", "CA", "C", "phi"}, // PHI: C0-N1-CA1-C1
  { 1, "N", "CA", "C", "N", "psi"}  // PSI: N0-CA0-C0-N1
};

void DihedralSearch::ListKnownTypes() {
  for (int i = 0; i < (int)NDIHTYPE; ++i)
    mprintf(" %s", DIH[i].name);
  mprintf("\n");
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
                                             const char* an0, const char* an1,
                                             const char* an2, const char* an3,
                                             const char* name) :
  offset_(off), aname0_(an0), aname1_(an1), aname2_(an2), aname3_(an3),
  name_(name)
{}

// DihedralSearch::DihedralToken::FindDihedralAtoms()
DihedralSearch::DihedralMask 
  DihedralSearch::DihedralToken::FindDihedralAtoms(Topology const& topIn, int resIn)
{
  int firstAtomRes = resIn;
  int lastAtomRes = resIn;
  if (offset_ == -1)
    --firstAtomRes;
  else if (offset_ == 1)
    ++lastAtomRes;
  int atom1 = topIn.FindAtomInResidue(firstAtomRes, aname0_);
  if (atom1 == -1) return DihedralMask();
  int atom2 = topIn.FindAtomInResidue(resIn, aname1_);
  if (atom2 == -1) return DihedralMask();
  int atom3 = topIn.FindAtomInResidue(resIn, aname2_);
  if (atom3 == -1) return DihedralMask();
  int atom4 = topIn.FindAtomInResidue(lastAtomRes, aname3_);
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
  mprintf("\tFound %u dihedrals.\n", dihedrals_.size());
  return 0;
}

// DihedralSearch::Clear()
void DihedralSearch::Clear() {
  dihedralTokens_.clear();
  dihedrals_.clear();
}

void DihedralSearch::PrintTypes() {
  for (std::vector<DihedralToken>::iterator dih = dihedralTokens_.begin();
                                            dih != dihedralTokens_.end(); ++dih)
  {
    mprintf(" %s", (*dih).Name().c_str());
  }
}
