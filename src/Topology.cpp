#include <stack> // For ParseMask
#include <algorithm> // sort
#include <sstream> // ostringstream
#include "Topology.h"
#include "CpptrajStdio.h"
// DEBUG
//#include <cmath> // sqrt

// CONSTRUCTOR
Topology::Topology() :
  debug_(0),
  topology_error_(0),
  firstSolventMol_(-1),
  NsolventMolecules_(0),
  finalSoluteRes_(-1),
  pindex_(0),
  nframes_(0),
  ntypes_(0),
  massptr_(0)
{ }

// DESTRUCTOR
Topology::~Topology() {
  if (massptr_!=0) delete[] massptr_;
}

// Topology::SetDebug()
void Topology::SetDebug(int debugIn) {
  debug_ = debugIn;
}

// Topology::SetParmName()
void Topology::SetParmName(std::string& title, const char* filename) {
  // NOTE: Check for NULL?
  parmName_ = title;
  fileName_.assign( filename );
}

// Topology::SetGBradiiSet()
void Topology::SetGBradiiSet(std::string &gbset) {
  radius_set_ = gbset;
}  

// Topology::SetPindex()
void Topology::SetPindex(int pindexIn) {
  pindex_ = pindexIn;
}

// Topology::SetReferenceCoords()
void Topology::SetReferenceCoords( Frame *frameptr ) {
    if (frameptr==NULL) return;
    refCoords_ = *frameptr;
}

// -----------------------------------------------------------------------------
// Topology::FinalSoluteRes()
/** Return 1 past the last solute residue. */
int Topology::FinalSoluteRes() {
  return finalSoluteRes_ + 1;
}

// Topology::Pindex()
int Topology::Pindex() {
  return pindex_;
}

// Topology::Natom()
int Topology::Natom() {
  return (int)atoms_.size();
}

// Topology::Nres()
int Topology::Nres() {
  return (int)residues_.size();
}

// Topology::Nmol()
int Topology::Nmol() {
  return (int)molecules_.size();
}

// Topology::FirstSolventMol()
int Topology::FirstSolventMol() {
  return firstSolventMol_;
}

// Topology::Nsolvent()
int Topology::Nsolvent() {
  return NsolventMolecules_;
}

// Topology::Nframes()
int Topology::Nframes() {
  return nframes_;
}

// Topology::Ntypes()
int Topology::Ntypes() {
  return ntypes_;
}

// Topology::IncreaseFrames()
void Topology::IncreaseFrames(int frames) {
  nframes_ += frames;
}

// Topology::c_str()
/** Return a printf-compatible char* of the parm name, or the full parm
  * filename if the parm name is empty.
  */
const char *Topology::c_str() {
  if (!parmName_.empty())
    return parmName_.c_str();
  return fileName_.c_str();
}

// Topology::ParmName()
std::string Topology::ParmName() { 
  return parmName_;
}

// Topology::OriginalFilename()
std::string Topology::OriginalFilename() {
  return fileName_;
}

// Topology::GBradiiSet()
std::string Topology::GBradiiSet() {
  return radius_set_;
}

// -----------------------------------------------------------------------------
// Topology::ResAtomStart()
Topology::atom_iterator Topology::ResAtomStart(int resnum) const {
  if (resnum < 0 || resnum >= (int)residues_.size())
    return atoms_.end();
  return atoms_.begin() + residues_[resnum].FirstAtom();
}

// Topology::ResAtomEnd()
Topology::atom_iterator Topology::ResAtomEnd(int resnum) const {
  if (resnum < 0 || resnum >= (int)residues_.size()-1)
    return atoms_.end();
  return atoms_.begin() + residues_[resnum+1].FirstAtom();
}

// Topology::MolAtomStart()
Topology::atom_iterator Topology::MolAtomStart(int molnum) const {
  if (molnum < 0 || molnum >= (int)molecules_.size())
    return atoms_.end();
  return atoms_.begin() + molecules_[molnum].BeginAtom();
}

// Topology::MolAtomEnd()
Topology::atom_iterator Topology::MolAtomEnd(int molnum) const {
  if (molnum < 0 || molnum >= (int) molecules_.size()-1)
    return atoms_.end();
  return atoms_.begin() + molecules_[molnum+1].BeginAtom();
}

// Topology::operator[]
const Atom& Topology::operator[](int idx) {
  return atoms_[idx];
}

// -----------------------------------------------------------------------------
// Topology::Res()
const Residue& Topology::Res(int idx) {
  return residues_[idx];
}

// -----------------------------------------------------------------------------
// Topology::SolventStart()
Topology::mol_iterator Topology::SolventStart() const {
  if (NsolventMolecules_==0)
    return molecules_.end();
  return molecules_.begin() + firstSolventMol_;
}

// Topology::SolventEnd()
Topology::mol_iterator Topology::SolventEnd() const {
  if (NsolventMolecules_==0)
    return molecules_.end();
  return molecules_.begin() + firstSolventMol_ + NsolventMolecules_;
}
// -----------------------------------------------------------------------------
int Topology::GetBondParamIdx( int idx, double &Rk, double &Req) {
  if (idx < 0 || idx > (int)bondrk_.size()) return 1;
  Rk = bondrk_[idx];
  Req = bondreq_[idx];
  return 0;
}
double Topology::GetBondedCutoff(int atom1, int atom2) {
  if (atom1 < 0 || atom2 < 0) return -1;
  if (atom1 >= Natom() || atom2 >= Natom()) return -1;
  return GetBondedCut( atoms_[atom1].Element(), atoms_[atom2].Element() );
}
// -----------------------------------------------------------------------------
// Topology::ResAtomRange()
int Topology::ResAtomRange(int res, int *resstart, int *resstop) {
  int nres1 = (int)residues_.size() - 1;
  if (res < 0 || res > nres1) return 1;
  *resstart = residues_[res].FirstAtom();
  if (res == nres1)
    *resstop = (int)atoms_.size();
  else
    *resstop = residues_[res+1].FirstAtom();
  return 0;
}

// Topology::ResidueName()
char *Topology::ResidueName(int res) {
  if (res < 0 || res >= (int) residues_.size())
    return NULL;
  return (char*)residues_[res].c_str();
}

// Topology::ResAtomName()
/** Given an atom number, set buffer with residue name and number along with
  * atom name with format: <resname[res]><res+1>@<atomname>, e.g. ARG_11@CA.
  * Replace any blanks in resname with '_'.
  */
std::string Topology::ResAtomName(int atom) {
  std::string res_name;
  if (atom < 0 || atom >= (int)atoms_.size()) return res_name;
  std::string atom_name( atoms_[atom].c_str() );
  int res = atoms_[atom].ResNum();
  res_name.assign( residues_[res].c_str() );
  // NOTE: ensure a residue size of 4?
  if (res_name[3]==' ')
    res_name[3]='_';
  ++res; // want output as res+1
  std::ostringstream oss;
  oss << res_name << res << "@" << atom_name;
  return oss.str();
}

// Topology::FindAtomInResidue()
/** Find the atom # of the specified atom name in the given residue.
  * \param res Residue number to search.
  * \param atname Atom name to find.
  * \return the atom number of the specified atom if found in the given residue.
  * \return -1 if atom not found in given residue.
  */
int Topology::FindAtomInResidue(int res, NameType atname) {
  if (res < 0 || res >= (int)residues_.size()) return -1;
  for (atom_iterator atom = ResAtomStart(res); atom != ResAtomEnd(res); atom++)
    if ( (*atom).Name() == atname )
      return ( atom - atoms_.begin() );
  return -1;
}

// Topology::FindResidueMaxNatom()
/** Return the # atoms in the largest residue. */
int Topology::FindResidueMaxNatom() {
  if (residues_.size() <= 1)
    return (int)atoms_.size();
  int largest_natom = 0;
  int lastatom = (int)atoms_.size();
  for (std::vector<Residue>::iterator res = residues_.end() - 1;
                                      res != residues_.begin(); res--)
  {
    int firstatom = (*res).FirstAtom();
    int diff = lastatom - firstatom;
    if (diff > largest_natom) largest_natom = diff;
    lastatom = firstatom;
  }
  return largest_natom;
}

// Topology::SoluteAtoms()
int Topology::SoluteAtoms() {
  if (NsolventMolecules_ == 0)
    return (int)atoms_.size();
  return molecules_[firstSolventMol_].BeginAtom();
}

// Topology::Mass() 
// NOTE: Stopgap - need to figure out a better way
// TODO: Figure out a better way to set up frames
double *Topology::Mass() {
  if (atoms_.empty()) return 0;
  if (massptr_ == 0) {
    massptr_ = new double[ atoms_.size() ];
    unsigned int m = 0;
    for (std::vector<Atom>::iterator atom = atoms_.begin();
                                     atom != atoms_.end(); atom++) 
      massptr_[m++] = (*atom).Mass();
  }
  return massptr_;
}
// -----------------------------------------------------------------------------
// Topology::Summary()
void Topology::Summary() {
  mprintf("\t\tTopology contains %zu atoms.\n", atoms_.size());
  mprintf("\t\t                  %zu residues.\n", residues_.size());
  mprintf("\t\t                  %zu bonds.\n", (bonds_.size()+bondsh_.size()) / 3 );
  mprintf("\t\t                  %zu molecules.\n", molecules_.size());
  mprintf("\t\t                  Box: %s\n",box_.TypeName());
  if (NsolventMolecules_>0) {
    mprintf("\t\t                  %i solvent molecules.\n", NsolventMolecules_);
    mprintf("\t\t                  First solvent mol is %i\n", firstSolventMol_+1);
    mprintf("\t\t                  Final solute residue is %i\n", finalSoluteRes_+1);
  }
  /*if (!bondsh_.empty() || !bonds_.empty())
    mprintf("  %zu bonds to hydrogen, %zu other bonds.\n",bondsh_.size()/3,
            bonds_.size()/3);*/
}

// Topology::ParmInfo()
void Topology::ParmInfo() {
  mprintf(" %i: %s, %zu atoms, %zu res, box: %s, %zu mol",pindex_,c_str(),
          atoms_.size(),residues_.size(),box_.TypeName(),molecules_.size());
  if (NsolventMolecules_>0)
    mprintf(", %i solvent",NsolventMolecules_);
  if (nframes_>0) 
    mprintf(", %i frames",nframes_);
  mprintf("\n");
}

// Topology::PrintAtomInfo()
void Topology::PrintAtomInfo(const char *maskString) {
  AtomMask mask;
  mask.SetMaskString( (char*)maskString ); // TODO: Use only strings
  ParseMask(refCoords_, mask, true);
  if ( mask.None() )
    mprintf("\tSelection is empty.\n");
  else
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); atom++) {
      mprintf("\tAtom %i", *atom+1);
      atoms_[ (*atom) ].Info();
    }
} 

// Topology::PrintBonds()
void Topology::PrintBonds(std::vector<int>& barray) {
  for (std::vector<int>::iterator batom = barray.begin();
                                  batom != barray.end(); batom++)
  {
    int atom1 = ((*batom) / 3);
    ++batom;
    int atom2 = ((*batom) / 3);
    mprintf("\tAtom %i:%s to %i:%s", atom1+1, atoms_[atom1].c_str(),
                                     atom2+1, atoms_[atom2].c_str());
    ++batom;
    if (*batom==-1) {
      double req = GetBondedCut(atoms_[atom1].Element(),atoms_[atom2].Element());
      mprintf("  EQ=%lf\n", req);
    } else {
      // TODO: Bond index should be -1
      double req = bondreq_[*batom - 1];
      double rk = bondrk_[*batom - 1];
      mprintf("  EQ=%lf K=%lf\n", req, rk);
    }
  }
}

// Topology::PrintBondInfo()
void Topology::PrintBondInfo() {
  if (!bondsh_.empty()) {
    mprintf("%zu BONDS TO HYDROGEN:\n",bondsh_.size()/3);
    PrintBonds( bondsh_ );
  }
  if (!bonds_.empty()) {
    mprintf("%zu BONDS TO NON-HYDROGEN:\n",bonds_.size()/3);
    PrintBonds( bonds_ );
  }
}

// Topology::PrintMoleculeInfo()
void Topology::PrintMoleculeInfo() {
  if (molecules_.empty())
    mprintf("\t[%s] No molecule info.\n",c_str());
  else {
    mprintf("MOLECULES:\n");
    unsigned int mnum = 1;
    for (std::vector<Molecule>::iterator mol = molecules_.begin(); 
                                         mol != molecules_.end(); mol++)
    {
      int firstres = (*mol).FirstRes();
      mprintf("\tMolecule %u, %i atoms, first residue %i:%s\n",mnum,(*mol).NumAtoms(),
              firstres+1,residues_[firstres].c_str());
      ++mnum;
    }
  }
}

// Topology::PrintResidueInfo()
void Topology::PrintResidueInfo() {
  unsigned int rnum = 1;
  for (std::vector<Residue>::iterator res = residues_.begin();
                                      res != residues_.end(); res++)
  {
    mprintf("\tResidue %u, %s, first atom %i\n",rnum,(*res).c_str(),(*res).FirstAtom()+1);
    ++rnum;
  }
}

// -----------------------------------------------------------------------------
// Topology::AddAtom()
void Topology::AddAtom(Atom atomIn, Residue resIn, const double* XYZin) {
  // DEBUG
  //mprintf("Adding atom %i [%s] of residue %i [%s]\n",
  //       atomIn.Num(), atomIn.Name(), resIn.Num(), resIn.Name());
  // Check if this is a new residue
  if (residues_.empty() || resIn != residues_.back()) {
    // Set atom number of first atom of this residue
    resIn.SetFirstAtom( atoms_.size() );
    residues_.push_back( resIn );
  }
  //mprintf("DEBUG: Atom %zu belongs to residue %zu\n",atoms_.size(), residues_.size());
  //mprintf("DEBUG: Atom %zu: %lf %lf %lf\n",atoms_.size(),XYZin[0],XYZin[1],XYZin[2]);
  // Overwrite atom number
  //atomIn.SetNum( atoms_.size() );
  atomIn.SetResNum( residues_.size()-1 );
  atoms_.push_back(atomIn);
  // Add coordinate if given
  refCoords_.AddXYZ( XYZin );
}

// Topology::StartNewMol()
void Topology::StartNewMol() {
  // If this is the first time this routine has been called, consider all
  // atoms to this point as belonging to first molecule. 
  if (molecules_.empty()) {
    //mprintf("DEBUG:\tFirst molecule, atoms 0 to %zu\n",atoms_.size());
    molecules_.push_back( Molecule(0, 0, atoms_.size()) );
  } else {
    // The first atom of this molecule will be end atom of last molecule.
    int molFirstAtom = molecules_.back().EndAtom();
    // Only add a new molecule if #atoms > first atom of the molecule.
    if ((int)atoms_.size() > molFirstAtom) {
      // Figure out first residue
      int molFirstRes = atoms_[molFirstAtom].ResNum();
      molecules_.push_back( Molecule(molFirstRes, molFirstAtom, atoms_.size()) );
    }
    // First atom
    //mprintf("DEBUG:\tMolecule %zu, atoms %i to %zu\n",
    //       molecules_.size(), lastAtom, atoms_.size());
  }
}

// Topology::CreateAtomArray()
int Topology::CreateAtomArray(std::vector<NameType>& names, std::vector<double>& charge,
                        std::vector<int>& at_num, std::vector<double>& mass,
                        std::vector<int>& atype_index, std::vector<NameType>& types,
                        std::vector<double>& gb_radii, std::vector<double>& gb_screen, 
                        std::vector<NameType>& resnames, std::vector<int>& resnums)
{
  if (names.empty() || resnames.empty()) {
    mprinterr("Error: Topology: Cannot create Atom/Residue arrays from empty input arrays.\n");
    return 1;
  }
  size_t natom = names.size();
  if ( natom != charge.size() ||
       natom != mass.size() ||
       natom != atype_index.size() ||
       natom != types.size()
     )
  {
    mprinterr("Error: Topology: Array sizes for #atoms from input parm do not match.\n");
    mprinterr("\tNames = %zu\n",names.size());
    mprinterr("\tCharges = %zu\n",charge.size());
    mprinterr("\tMass = %zu\n",mass.size());
    mprinterr("\tAtomType Index = %zu\n",atype_index.size());
    mprinterr("\tTypes = %zu\n",types.size());
    return 1;
  }
  if (resnames.size() != resnums.size()) {
    mprinterr("Error: Topology: Array sizes for #residues from input parm do not match.\n");
    mprinterr("\tResNames = %zu\n",resnames.size());
    mprinterr("\tResNums = %zu\n",resnums.size());
    return 1;
  }
  // ATOMIC_NUMBER may not be present
  if (at_num.empty()) {
    mprintf("Warning: [%s] ATOMIC_NUMBER not present in topology.\n", c_str());
    at_num.resize(natom, 0);
  }
  // GB params may be empty in old amber parm
  if (gb_radii.empty()) {
    mprintf("Warning: [%s] GB RADII not present in topology.\n", c_str());
    gb_radii.resize(natom, 0);
  }
  if (gb_screen.empty()) {
    mprintf("Warning: [%s] GB SCREEN not present in topology.\n", c_str());
    gb_screen.resize(natom, 0);
  }
  // Create atom information
  atoms_.reserve( natom );
  int resnum = 0;
  std::vector<int>::iterator Res = resnums.begin() + 1;
  for (size_t atom = 0; atom < natom; atom++) {
    if (Res!=resnums.end() && atom >= (size_t)*Res) {
      ++resnum;
      ++Res;
    }
    atoms_.push_back( Atom( names[atom], charge[atom], at_num[atom],
                            mass[atom], atype_index[atom], types[atom],
                            gb_radii[atom], gb_screen[atom], resnum)
                    );
  }
  // Create residue information
  size_t nres = resnames.size();
  residues_.reserve( nres );
  for (size_t res = 0; res < nres; res++)
    residues_.push_back( Residue( resnames[res], resnums[res] ) ); 

  return 0;
}

// Topology::CreateMoleculeArray()
int Topology::CreateMoleculeArray(std::vector<int> &atomsPerMol, Box parmbox, 
                                  int finalSoluteRes, int firstSolvMol)
{
  if (atomsPerMol.empty() || finalSoluteRes < 0 || firstSolvMol < 0) {
    mprinterr("Error: Topology: Could not set up molecule info.\n");
    mprinterr("\tAtomsPerMol size = %zu\n",atomsPerMol.size());
    mprinterr("\tfinalSoluteRes = %i\n",finalSoluteRes);
    mprinterr("\tfirstSolvMol = %i\n",firstSolvMol);
    return 1;
  }
  if (atoms_.empty()) {
    mprinterr("Error: Topology: Cannot set up molecule info; no atoms present.\n");
    return 1;
  }
  molecules_.clear();
  molecules_.reserve( atomsPerMol.size() );
  int molbegin = 0;
  int molend = 0;
  int molnum = 0;
  for (std::vector<int>::iterator molsize = atomsPerMol.begin(); 
                                  molsize != atomsPerMol.end(); molsize++)
  {
    int firstRes = atoms_[molbegin].ResNum();
    molend += *molsize;
    // Update atoms molecule numbers
    for (int at = molbegin; at < molend; at++)
      atoms_[at].SetMol( molnum );
    molecules_.push_back( Molecule(firstRes, molbegin, molend) );
    molbegin = molend;
    ++molnum;
  }
  box_ = parmbox;
  finalSoluteRes_ = finalSoluteRes;
  firstSolventMol_ = firstSolvMol;
  return 0;
}

// Topology::SetBondInfo()
int Topology::SetBondInfo(std::vector<int> &bonds, std::vector<int> &bondsh,
                          std::vector<double> &bond_rk, std::vector<double> &bond_req) 
{
  if (bonds.empty() && bondsh.empty()) {
    mprinterr("Error: Topology: Input bonds and bondsh are empty.\n");
    return 1;
  }
  if (atoms_.empty()) {
    mprinterr("Error: Topology: Cannot set up bonds, no atoms present.\n");
    return 1;
  }
  // NOTE: Clear any previous bond info here?
  bonds_ = bonds;
  bondsh_ = bondsh;
  SetAtomBondInfo();
  // Create bond parameter arrays
  if (bond_rk.size() != bond_req.size()) {
    mprinterr("Error: Topology: Bond parameters have different lengths (%zu != %zu)\n",
              bond_rk.size(), bond_req.size());
    return 1;
  }
  bondrk_ = bond_rk;
  bondreq_ = bond_req;

  /*else {
    std::vector<double>::iterator req = bond_req.begin();
    for (std::vector<double>::iterator rk = bond_rk.begin(); rk != bond_rk.end(); rk++)
    {
      bondParm_.push_back( ParmBondType( *rk, *req) );
      ++req;
    }
  }*/
  return 0;
}

// Topology::SetAngleInfo()
int Topology::SetAngleInfo(std::vector<int>& angles, std::vector<int>& anglesh,
                           std::vector<double>& tk, std::vector<double>& teq)
{
  angles_ = angles;
  anglesh_ = anglesh;
  angletk_ = tk;
  angleteq_ = teq;
  return 0;
}

// Topology::SetDihedralInfo()
int Topology::SetDihedralInfo(std::vector<int>& dihedrals, std::vector<int>& dihedralsh,
                              std::vector<double>& pk, std::vector<double>& pn,
                              std::vector<double>& phase,
                              std::vector<double>& scee, std::vector<double>& scnb)
{
  dihedrals_ = dihedrals;
  dihedralsh_ = dihedralsh;
  dihedralpk_ = pk;
  dihedralpn_ = pn;
  dihedralphase_ = phase;
  scee_ = scee;
  scnb_ = scnb;
  if (scee_.empty())
    scee_.resize( dihedralpk_.size(), 1.2);
  if (scnb.empty())
    scnb_.resize( dihedralpk_.size(), 2.0);
  return 0;
}

// Topology::SetAmberHbond()
int Topology::SetAmberHbond(std::vector<double>& asol, std::vector<double>& bsol,
                            std::vector<double>& hbcut)
{
  asol_ = asol;
  bsol_ = bsol;
  hbcut_ = hbcut;
  return 0;
}

// Topology::SetAmberExtra()
// TODO: Auto generate
int Topology::SetAmberExtra(std::vector<double>& solty, std::vector<NameType>& itree, 
                            std::vector<int>& join, std::vector<int>& irotat)
{
  solty_ = solty;
  itree_ = itree;
  join_ = join;
  irotat_ = irotat;
  return 0;
}

// Topology::SetNonbondInfo()
int Topology::SetNonbondInfo(int ntypesIn, std::vector<int>& nbindex, 
                             std::vector<double>& lja, std::vector<double>& ljb)
{
  ntypes_ = ntypesIn;
  int nbsize = ntypes_ * ntypes_;
  if ((int)nbindex.size() != nbsize) {
    mprinterr("Error: Topology: Size of NB_index (%zu) does not match ntypes*ntypes (%i)\n",
              nbindex.size(), nbsize);
    return 1;
  }
  int ljsize = (ntypes_ * (ntypes_+1)) / 2;
  if (ljsize != (int)lja.size() ||
      ljsize != (int)ljb.size()) 
  {
    mprinterr("Error: Topology: LJ parameters have wrong size, %i != (%zu | %zu)\n",
              ljsize, lja.size(), ljb.size());
    return 1;
  }
  nbindex_ = nbindex;
  lja_ = lja;
  ljb_ = ljb;
/*  std::vector<double>::iterator A = lja.begin();
  for (std::vector<double>::iterator B = ljb.begin(); B != ljb.end(); ljb++) {
    nonbondParm_.push_back( ParmNonbondType( *A, *B ) );
    ++A;
  }*/
  return 0;
}

// Topology::SetAtomBondInfo()
/** Set up bond information in the atoms array based on bonds and bondsh.
  */
void Topology::SetAtomBondInfo() {
  for (std::vector<int>::iterator bond = bonds_.begin(); bond != bonds_.end(); bond++) {
    int atom1 = *bond / 3;
    ++bond;
    int atom2 = *bond / 3;
    ++bond;
    atoms_[atom1].AddBond( atom2 );
    atoms_[atom2].AddBond( atom1 );
  }

  for (std::vector<int>::iterator bond = bondsh_.begin(); bond != bondsh_.end(); bond++) {
    int atom1 = *bond / 3;
    ++bond;
    int atom2 = *bond / 3;
    ++bond;
    atoms_[atom1].AddBond( atom2 );
    atoms_[atom2].AddBond( atom1 );
  }
}

// Topology::CommonSetup()
void Topology::CommonSetup(bool bondsearch, bool molsearch) {
  // Replace asterisks in atom names with single quote
  for (std::vector<Atom>::iterator atom = atoms_.begin(); atom != atoms_.end(); atom++)
    (*atom).ReplaceAsterisk();
  // Add placeholder to residues to indicate last residue
  //residues_.push_back( Residue( atoms_.size() ) );
  // Set up bond information if specified and necessary
  if (bondsearch) {
    if (bonds_.empty() && bondsh_.empty() && !refCoords_.empty()) {
      GetBondsFromAtomCoords();
      DetermineMolecules();
    }
  }
  
  if (molecules_.empty()) 
    DetermineMolecules();

  // Set up solvent information
  SetSolventInfo();

  // Determine excluded atoms
  DetermineExcludedAtoms(); 

}

// -----------------------------------------------------------------------------
// Topology::compareElement()
/** Compare a pair of elements a1 and a2 with target values b1 and b2.
  * If either combination of a1 and a2 match b1 and b2 (i.e. a1==b1 and 
  * a2==b2, or a1==b2 and a2==b1) return true, otherwise return false.
  */
bool Topology::compareElement(Atom::AtomicElementType a1, Atom::AtomicElementType a2,
                              Atom::AtomicElementType b1, Atom::AtomicElementType b2)
{
  if      (a1==b1 && a2==b2)
    return true;
  else if (a1==b2 && a2==b1)
    return true;
  return false;
}

// Topology::GetBondedCut() 
/** Return a cutoff based on optimal covalent bond distance based on the 
  * identities of atom1 and atom2. When multiple hybridizations are possible
  * the longest possible bond length is used.
  * Unless otherwise noted values taken from:
  * - Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 
  *       2nd ed., Butterworths, London, 1958; 
  * - B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, 
  *       No. 31, Washington, DC, 1970; S.W. Benson, J. Chem. Educ., 42, 502 (1965).
  * Can be found on the web at:
  * - http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
  */
// NOTE: Store cutoff^2 instead??
double Topology::GetBondedCut(Atom::AtomicElementType atom1, Atom::AtomicElementType atom2) 
{
  // Default cutoff
  double cut = 1.60;
  // Self
  if (atom1==atom2) {
    if      (atom1==Atom::HYDROGEN   ) cut=0.74;
    else if (atom1==Atom::NITROGEN   ) cut=1.45;
    else if (atom1==Atom::CARBON     ) cut=1.54;
    else if (atom1==Atom::OXYGEN     ) cut=1.48;
    else if (atom1==Atom::PHOSPHORUS ) cut=2.21;
    else if (atom1==Atom::SULFUR     ) cut=2.05; // S-S gas-phase value; S=S is 1.49
  }
  // Bonds to H 
  else if ( compareElement(atom1,atom2,Atom::HYDROGEN,Atom::CARBON) )
    cut=1.09;
  else if ( compareElement(atom1,atom2,Atom::HYDROGEN,Atom::NITROGEN) )
    cut=1.01;
  else if ( compareElement(atom1,atom2,Atom::HYDROGEN,Atom::OXYGEN) )
    cut=0.96;
  else if ( compareElement(atom1,atom2,Atom::HYDROGEN,Atom::PHOSPHORUS) )
    cut=1.44;
  else if ( compareElement(atom1,atom2,Atom::HYDROGEN,Atom::SULFUR) )
    cut=1.34;
  // Bonds to C
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::NITROGEN) )
    cut=1.47;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::OXYGEN) )
    cut=1.43;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::PHOSPHORUS) )
    cut=1.84;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::FLUORINE) )
    cut=1.35;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::CHLORINE) )
    cut=1.77;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::BROMINE) )
    cut=1.94;
  else if ( compareElement(atom1,atom2,Atom::CARBON,Atom::SULFUR) )
    cut=1.82;
  // Bonds to N
  else if ( compareElement(atom1,atom2,Atom::NITROGEN,Atom::OXYGEN) )
    cut=1.40;
  else if ( compareElement(atom1,atom2,Atom::NITROGEN,Atom::SULFUR) )
    cut=1.68; // Postma & Vos, Acta Cryst. (1973) B29, 915
  else if ( compareElement(atom1,atom2,Atom::NITROGEN,Atom::FLUORINE) )
    cut=1.36;
  else if ( compareElement(atom1,atom2,Atom::NITROGEN,Atom::CHLORINE) )
    cut=1.75;
  else if ( compareElement(atom1,atom2,Atom::NITROGEN,Atom::PHOSPHORUS) )
    cut=1.71; // Avg over all nX-pX from gaff.dat
  // Bonds to P
  else if ( compareElement(atom1,atom2,Atom::PHOSPHORUS,Atom::OXYGEN) )
    cut=1.63;
  else if ( compareElement(atom1,atom2,Atom::PHOSPHORUS,Atom::SULFUR) )
    cut=1.86;
  else if ( compareElement(atom1,atom2,Atom::PHOSPHORUS,Atom::FLUORINE) )
    cut=1.54;
  else if ( compareElement(atom1,atom2,Atom::PHOSPHORUS,Atom::CHLORINE) )
    cut=2.03;
  // Bonds to O
  else if ( compareElement(atom1,atom2,Atom::OXYGEN,Atom::SULFUR) )
    cut=1.48;
  else if ( compareElement(atom1,atom2,Atom::OXYGEN,Atom::FLUORINE) )
    cut=1.42;
  // Bonds to S
  else if ( compareElement(atom1,atom2,Atom::SULFUR,Atom::FLUORINE) )
    cut=1.56;
  else if ( compareElement(atom1,atom2,Atom::SULFUR,Atom::CHLORINE) )
    cut=2.07;
  // No cutoff, use default
  else {
    mprintf("Warning: GetBondedCut: Cut not found for %s - %s\n",
            Atom::AtomicElementName[atom1],Atom::AtomicElementName[atom2]);
    mprintf("                       Using default cutoff of %lf\n",cut);
  }

  //mprintf("\t\tCUTOFF: [%s] -- [%s] = %lf\n",AtomicElementName[atom1],
  //        AtomicElementName[atom2],cut);

  // Padding value
  cut += 0.15;
  return cut;
}

// Topology::GetBondsFromAtomCoords()
void Topology::GetBondsFromAtomCoords() {
  int stopatom;
  // Dont want to check for bonds between residues once we are in a 
  // solvent region, so set the last residue to be searched to the final
  // solute residue.
  std::vector<Residue>::iterator firstSolventResidue = residues_.end();

  mprintf("\t%s: determining bond info from distances.\n",c_str());
  // ----- STEP 1: Determine bonds within residues
  //int resnum = 0; // DEBUG 
  for (std::vector<Residue>::iterator res = residues_.begin(); 
                                      res != residues_.end(); res++) 
  {
    // Check if this residue is the first solvent residue.
    if (firstSolventResidue!=residues_.end() && (*res).NameIsSolvent())
      firstSolventResidue = res;
    // Get residue start atom.
    int startatom = (*res).FirstAtom();
    // Get residue end atom.
    std::vector<Residue>::iterator nextres = res + 1;
    if (nextres == residues_.end())
      stopatom = atoms_.size();
    else
      stopatom = (*nextres).FirstAtom();
    // DEBUG
    //mprintf("\tRes %i Start atom %zu coords: ",resnum+1, startatom+1);
    //atoms_[startatom].PrintXYZ();
    //mprintf("\n");
    // Check for bonds between each atom in the residue.
    for (int atom1 = startatom; atom1 < stopatom - 1; atom1++) {
      // If this is a hydrogen and it already has a bond, move on.
      if (atoms_[atom1].Element()==Atom::HYDROGEN &&
          atoms_[atom1].Nbonds() > 0 )
        continue;
      for (int atom2 = atom1 + 1; atom2 < stopatom; atom2++) {
        double D2 = refCoords_.DIST2( atom1, atom2 );
        double cutoff2 = GetBondedCut(atoms_[atom1].Element(), atoms_[atom2].Element());
        //mprintf("\t\t%i:[%s] -- %i:[%s] D=%lf  Cut=%lf\n",atom1+1,atoms_[atom1].c_str(),
        //         atom2+1,atoms_[atom2].c_str(),sqrt(D2),cutoff2);
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) {
          AddBond(atom1, atom2);
          // Once a bond has been made to hydrogen move on.
          if (atoms_[atom1].Element()==Atom::HYDROGEN) break;
        }
      }
    }
    //++resnum; // DEBUG
  }

  // ----- STEP 2: Determine bonds between adjacent residues
  std::vector<Molecule>::iterator nextmol = molecules_.begin();
  if (!molecules_.empty())
    ++nextmol;
  for (std::vector<Residue>::iterator res = residues_.begin() + 1;
                                      res != firstSolventResidue; res++)
  {
    // If molecule information is already present, check if first atom of 
    // this residue >= first atom of next molecule, which indicates this
    // residue and the previous residue are in different molecules.
    if ( (nextmol != molecules_.end()) && 
         ((*res).FirstAtom() >= (*nextmol).BeginAtom()) )
    {
      ++nextmol;
      continue;
    }
    // Get first residue start atom
    std::vector<Residue>::iterator previous_res = res - 1;
    int startatom = (*previous_res).FirstAtom();
    // First residue stop atom, second residue start atom
    int midatom = (*res).FirstAtom();
    // Second residue stop atom
    std::vector<Residue>::iterator nextres = res + 1;
    if (nextres == residues_.end())
      stopatom = atoms_.size();
    else
      stopatom = (*nextres).FirstAtom();
    //mprintf("\tBonds between residues %s and %s\n",(*previous_res).c_str(),(*res).c_str());
    // Check for bonds between adjacent residues
    for (int atom1 = startatom; atom1 < midatom; atom1++) {
      if (atoms_[atom1].Element()==Atom::HYDROGEN) continue;
      for (int atom2 = midatom; atom2 < stopatom; atom2++) {
        if (atoms_[atom2].Element()==Atom::HYDROGEN) continue;
        double D2 = refCoords_.DIST2( atom1, atom2 );
        double cutoff2 = GetBondedCut(atoms_[atom1].Element(), atoms_[atom2].Element());
        //mprintf("\t\t%i:[%s] -- %i:[%s] D=%lf  Cut=%lf\n",atom1+1,atoms_[atom1].c_str(),
        //         atom2+1,atoms_[atom2].c_str(),sqrt(D2),cutoff2);
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) 
          AddBond(atom1, atom2);
      }
    }
  }

  mprintf("\t%s: %zu bonds to hydrogen, %zu other bonds.\n",c_str(),
          bondsh_.size()/3, bonds_.size()/3);
}

// Topology::ClearBondInfo()
void Topology::ClearBondInfo() {
  bonds_.clear();
  bondsh_.clear();
  for (std::vector<Atom>::iterator atom = atoms_.begin(); atom != atoms_.end(); atom++)
    (*atom).ClearBonds();
}

// Topology::AddBond()
/** Create a bond between atom1 and atom2, update the atoms array.
  * For bonds to H always insert the H second.
  */
void Topology::AddBond(int atom1, int atom2) {
  bool a1H = (atoms_[atom1].Element() == Atom::HYDROGEN);
  bool a2H = (atoms_[atom2].Element() == Atom::HYDROGEN);
  //mprintf("\t\t\tAdding bond %i to %i (isH=%i)\n",atom1+1,atom2+1,(int)isH);
  // Update bonds arrays
  // TODO: Check for duplicates
  if (a1H || a2H) {
    if (a1H) {
      bondsh_.push_back(atom2*3);
      bondsh_.push_back(atom1*3);
    } else {
      bondsh_.push_back(atom1*3);
      bondsh_.push_back(atom2*3);
    }
    bondsh_.push_back(-1);
  } else  {
    bonds_.push_back(atom1*3);
    bonds_.push_back(atom2*3);
    bonds_.push_back(-1);
  }
  // Update atoms
  atoms_[atom1].AddBond( atom2 );
  atoms_[atom2].AddBond( atom1 );
}

// Topology::VisitAtom()
// NOTE: Use iterator instead of atom num?
void Topology::VisitAtom(int atomnum, int mol) {
  // Return if this atom already has a molecule number
  if (!atoms_[atomnum].NoMol()) return;
  // Mark this atom as visited
  atoms_[atomnum].SetMol( mol );
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = atoms_[atomnum].bondbegin();
                           bondedatom != atoms_[atomnum].bondend(); bondedatom++)
    VisitAtom(*bondedatom, mol);
}

// Topology::DetermineMolecules()
/** Determine individual molecules using bond information. Performs a 
  * recursive search over the bonds of each atom.
  */
void Topology::DetermineMolecules() {
  std::vector<Atom>::iterator atom;

  mprintf("\t%s: determining molecule info from bonds.\n",c_str());
  // Reset molecule info for each atom
  for (atom = atoms_.begin(); atom != atoms_.end(); atom++)
    (*atom).SetMol( -1 );
  // Perform recursive search along bonds of each atom
  int mol = 0;
  int atomnum = 0;
  for (atom = atoms_.begin(); atom != atoms_.end(); atom++)
  {
    if ( (*atom).NoMol() ) {
      VisitAtom( atomnum, mol );
      ++mol;
    }
    ++atomnum;
  }
  if (debug_>0)
    mprintf("\t%i molecules.\n",mol);

  // Update molecule information
  molecules_.resize( mol );
  if (mol == 0) return;
  std::vector<Molecule>::iterator molecule = molecules_.begin();
  (*molecule).SetFirst(0, 0);
  atom = atoms_.begin(); 
  int lastMol = (*atom).Mol();
  int atomNum = 0;
  for (; atom != atoms_.end(); atom++)
  {
    if ( (*atom).Mol() != lastMol ) {
      // Set last atom of molecule
      (*molecule).SetLast( atomNum );
      // Set first atom and resnum of next molecule
      ++molecule;
      (*molecule).SetFirst( atomNum, (*atom).ResNum() );
      lastMol = (*atom).Mol();
    }
    ++atomNum;
  }
  (*molecule).SetLast( atoms_.size() );
}

// Topology::AtomDistance()
void Topology::AtomDistance(int originalAtom, int atom, int dist, std::set<int> &excluded) 
{
  // If this atom is already too far away return
  if (dist==4) return;
  // dist is less than 4 and this atom greater than original, add exclusion
  if (atom > originalAtom)
    excluded.insert( atom ); 
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = atoms_[atom].bondbegin();
                           bondedatom != atoms_[atom].bondend();
                           bondedatom++)
    AtomDistance(originalAtom, *bondedatom, dist+1, excluded);
}

// Topology::DetermineExcludedAtoms()
/** For each atom, determine which atoms with greater atom# are within
  * 4 bonds (and therefore should be excluded from a non-bonded calc).
  */
void Topology::DetermineExcludedAtoms() {
  // A set is used since it automatically sorts itself and rejects duplicates.
  std::set<int> excluded_i;
  int natom = (int)atoms_.size();
  for (int atomi = 0; atomi < natom; atomi++) {
    excluded_i.clear();
    //mprintf("\tDetermining excluded atoms for atom %i\n",atomi+1);
    // AtomDistance recursively sets each atom bond distance from atomi
    AtomDistance(atomi, atomi, 0, excluded_i);
    atoms_[atomi].AddExclusionList( excluded_i );
    // DEBUG
    //mprintf("Atom %i Excluded:");
    //for (std::set<int>::iterator ei = excluded_i.begin(); ei != excluded_i.end(); ei++)
    //  mprintf(" %i",*ei + 1);
    //mprintf("\n");
  } // END loop over atomi
}

// -----------------------------------------------------------------------------
// Topology::SetSolventInfo()
/** Determine which molecules are solvent. Set finalSoluteRes and 
  * firstSolventMol if not already set.
  */
void Topology::SetSolventInfo() {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolventInfo: No molecule information.\n");
    topology_error_ = 1;
    return;
  }
  NsolventMolecules_ = 0;
  int numSolvAtoms = 0;

  if (firstSolventMol_ == -1 ) { // NOTE: Not checking firstSoluteRes
    // Loop over each molecule. Check if first residue of molecule
    // is solvent.
    finalSoluteRes_ = -1;
    for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                         mol != molecules_.end(); mol++)
    {
      int firstRes = (*mol).FirstRes();
      if ( residues_[firstRes].NameIsSolvent() ) {
        (*mol).SetSolvent();
        ++NsolventMolecules_;
        numSolvAtoms += (*mol).NumAtoms();
        if (firstSolventMol_==-1) {
          // Final solute residue is the one prior to this
          //finalSoluteRes = residues_.begin() + firstRes - 1;
          finalSoluteRes_ = firstRes - 1;
          // This is the first solvent molecule
          firstSolventMol_ = (int)(mol - molecules_.begin());
        }
      }
    }
  } else {
    // Every molecule from firstSolventMol on is solvent.
    for (std::vector<Molecule>::iterator mol = molecules_.begin() + firstSolventMol_;
                                         mol != molecules_.end(); mol++)
    {
      (*mol).SetSolvent();
      ++NsolventMolecules_;
      numSolvAtoms += (*mol).NumAtoms();
    }
  }
  if (firstSolventMol_==-1 && finalSoluteRes_ == -1)
    finalSoluteRes_ = (int)residues_.size() - 1;
  if (debug_>0) {
    if (NsolventMolecules_ == 0) 
      mprintf("    No solvent.\n");
    else
      mprintf("    %i solvent molecules, %i solvent atoms\n",NsolventMolecules_,numSolvAtoms);
  }
}

// -----------------------------------------------------------------------------
// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask) { 
  return ParseMask(refCoords_, mask, true);
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask) {
  return ParseMask(refCoords_, mask, false);
}

// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask, Frame &frame) {
  return ParseMask( frame, mask, true );
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask, Frame &frame) {
  return ParseMask( frame, mask, false );
}

// Topology::Mask_SelectDistance()
void Topology::Mask_SelectDistance( Frame &REF, char *mask, bool within, 
                                    bool byAtom, double distance ) 
{
  int endatom;
  if (REF.empty()) {
    mprinterr("Error: No reference set for [%s], cannot select by distance.\n",c_str());
    return;
  }
  // Distance has been pre-squared.
  // Create temporary array of atom #s currently selected in mask. Also
  // reset mask, it will be the output mask.
  std::vector<unsigned int> selected;
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask[i]=='T') {
      selected.push_back( i );
      mask[i] = 'F';
    }
  }
  if (selected.empty()) {
    mprinterr("Error: SelectAtomsWithin(%lf): No atoms in prior selection.\n",distance);
    return;
  }
/*  if (debug_ > 1) {
    mprintf("\t\t\tDistance Op: Within=%i  byAtom=%i  distance^2=%lf\n",
            (int)within, (int)byAtom, distance);
    mprintf("\t\t\tInitial Mask=[");
    for (std::vector<unsigned int>::iterator at = selected.begin(); at != selected.end(); at++)
      mprintf(" %u",*at + 1);
    mprintf(" ]\n");
  }*/

  if (byAtom) { // Select by atom
    // Loop over all atoms
    // TODO: OpenMP parallelize
    for (int atomi = 0; atomi < (int)atoms_.size(); atomi++) {
      // No need to calculate if atomi already selected
      if (mask[atomi] == 'T') continue;
      // Loop over initially selected atoms
      for (int idx = 0; idx < (int)selected.size(); idx++) {
        int atomj = selected[idx];
        double d2 = REF.DIST2(atomi, atomj);
        if (within) {
          if (d2 < distance) {
            mask[atomi] = 'T';
            break;
          }
        } else {
          if (d2 > distance) {
            mask[atomi] = 'T';
            break;
          }
        }
      }
    }
  } else { // Select by residue
    // TODO: OpenMP parallelize
    bool selectresidue = false;
    for (int atomi = 0; atomi < (int)atoms_.size(); atomi++) {
      // No need to calculate if atomi already selected
      if (mask[atomi] == 'T') continue;
      // Loop over initially selected atoms
      for (int idx = 0; idx < (int)selected.size(); idx++) {
        int atomj = selected[idx];
        //mprintf("\t\t\tAtom %i to atom %i\n",atomi+1,atomj+1);
        double d2 = REF.DIST2(atomi, atomj);
        if (within) {
          if (d2 < distance) selectresidue = true;
        } else {
          if (d2 > distance) selectresidue = true;
        }
        if (selectresidue) {
          int currentres = atoms_[atomi].ResNum();
          if (currentres+1 == (int)residues_.size())
            endatom = atoms_.size();
          else
            endatom = residues_[currentres+1].FirstAtom();
          //mprintf("\t\t\t\tSelecting residue %i (%i-%i)\n",currentres+1,
          //        residues_[currentres].FirstAtom()+1,endatom);
          for (int resatom = residues_[currentres].FirstAtom(); resatom < endatom; resatom++)
            mask[resatom] = 'T';
          selectresidue = false;
          break;
        }
      }
    }
  }

}

// Topology::Mask_AND()
void Topology::Mask_AND(char *mask1, char *mask2) {
  //mprintf("\t\t\tPerforming AND on masks.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    //mprintf(" [%c|%c]",mask1[i],mask2[i]);
    if (mask1[i]=='F' || mask2[i]=='F')
      mask1[i] = 'F';
    // Otherwise mask1 should already be T
  }
  //mprintf("\n");
}

// Topology::Mask_OR()
void Topology::Mask_OR(char *mask1, char *mask2) {
  //mprintf("\t\t\tPerforming OR on masks.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T' || mask2[i]=='T')
      mask1[i] = 'T';
    else
      mask1[i] = 'F';
  }
}

// Topology::Mask_NEG()
void Topology::Mask_NEG(char *mask1) {
  //mprintf("\t\t\tNegating mask.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T')
      mask1[i] = 'F';
    else
      mask1[i] = 'T';
  }
}

// Topology::MaskSelectResidues()
void Topology::MaskSelectResidues(NameType name, char *mask) {
  int endatom;
  std::vector<Residue>::iterator res1 = residues_.begin() + 1;

  //mprintf("\t\t\tSelecting residues named [%s]\n",*name);
  for (std::vector<Residue>::iterator res = residues_.begin();
                                      res != residues_.end(); res++)
  {
    if ( (*res).Name().Match( name ) ) {
      if (res1 == residues_.end())
        endatom = atoms_.size();
      else
        endatom = (*res1).FirstAtom();
      std::fill(mask + (*res).FirstAtom(), mask + endatom, 'T');
    }
    ++res1;
  }
}

// Topology::MaskSelectResidues()
// Mask args expected to start from 1
void Topology::MaskSelectResidues(int res1, int res2, char *mask) {
  int startatom, endatom;
  int nres = (int) residues_.size();
  //mprintf("\t\t\tSelecting residues %i to %i\n",res1,res2);
  // Get start atom
  if (res1 > nres) {
    if (debug_>0)
      mprintf("Warning: Select residues: res 1 out of range (%i)\n",res1);
    return;
  }
  startatom = residues_[res1-1].FirstAtom();
  // Get end atom
  //if ( res2 > nres) {
  //  mprinterr("Error: Select residues: res 2 out of range (%i)\n",res2);
  //  return;
  //}
  if ( res2 >= nres )
    endatom = (int)atoms_.size();
  else
    endatom = residues_[res2].FirstAtom();
  // Select atoms
  std::fill(mask + startatom, mask + endatom, 'T');
}

// Topology::MaskSelectAtoms()
void Topology::MaskSelectAtoms( NameType name, char *mask) {
  //mprintf("\t\t\tSelecting atoms named [%s]\n",*name);
  unsigned int m = 0;
  for (std::vector<Atom>::iterator atom = atoms_.begin();
                                   atom != atoms_.end(); atom++)
  {
    //mprintf("\t\t\t%u PARM[%s]  NAME[%s]",m,(*atom).c_str(),*name);
    if ( (*atom).Name().Match( name ) )
      mask[m] = 'T';
    //mprintf(" %c\n",mask[m]);
    ++m;
  } 
}

// Topology::MaskSelectAtoms()
// Mask args expected to start from 1
void Topology::MaskSelectAtoms(int atom1, int atom2, char *mask) {
  int startatom, endatom;
  //mprintf("\t\t\tSelecting atoms %i to %i\n",atom1,atom2);
  if (atom1 > (int)atoms_.size()) {
    if (debug_>0) 
      mprintf("Warning: Select atoms: atom 1 out of range (%i)\n",atom1);
    return;
  }
  startatom = atom1 - 1;
  if (atom2 > (int)atoms_.size()) 
    //mprinterr("Error: Select atoms: atom 2 out of range (%i)\n",atom2)
    endatom = atoms_.size();
  else
    endatom = atom2;
  // Select atoms
  std::fill(mask + startatom, mask + endatom, 'T');
}

// Topology::ParseMask()
bool Topology::ParseMask(Frame &REF, AtomMask &maskIn, bool intMask) {
  std::stack<char*> Stack;
  char *pMask = NULL; 
  char *pMask2 = NULL;

  for (AtomMask::token_iterator token = maskIn.begintoken();
                                token != maskIn.endtoken(); token++)
  {
    if (pMask==NULL) {
      // Create new blank mask
      pMask = new char[ atoms_.size() ];
      std::fill(pMask, pMask + atoms_.size(), 'F');
    }
    switch ( (*token).Type() ) {
      case MaskToken::ResNum : 
        MaskSelectResidues( (*token).Res1(), (*token).Res2(), pMask );
        break;
      case MaskToken::ResName :
        MaskSelectResidues( (*token).Name(), pMask );
        break;
      case MaskToken::AtomNum :
        MaskSelectAtoms( (*token).Res1(), (*token).Res2(), pMask );
        break;
      case MaskToken::AtomName :
        MaskSelectAtoms( (*token).Name(), pMask );
        break;
      case MaskToken::SelectAll :
        std::fill(pMask, pMask + atoms_.size(), 'T');
        break;
      case MaskToken::OP_AND :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_AND( Stack.top(), pMask2 );
        delete[] pMask2;
        break;
      case MaskToken::OP_OR :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_OR( Stack.top(), pMask2 );
        delete[] pMask2;
        break;
      case MaskToken::OP_NEG :
        Mask_NEG( Stack.top() );
        break;
      case MaskToken::OP_DIST :
        Mask_SelectDistance( REF, Stack.top(), (*token).Within(), (*token).ByAtom(), 
                             (*token).Distance() );
        break;
      default:
        mprinterr("Error: Invalid mask token (Mask [%s]).\n",maskIn.MaskString());
    }
    // Check if this mask should now go on the stack
    if ( (*token).OnStack() ) {
      //mprintf("Pushing Mask on stack, last Token [%s]\n",(*token).TypeName());
      Stack.push( pMask );
      pMask = NULL;
    }
  }
  // If pMask is not NULL it is probably a blank leftover
  if (pMask!=NULL) delete[] pMask;

  // If stack is empty here there was an error.
  if (Stack.empty()) {
    mprinterr("Error: Could not parse mask [%s].\n",maskIn.MaskString());
    return true;
  }

  // Top of the stack should point to the final mask
  pMask = Stack.top();
  Stack.pop();
  // Stack should be empty now
  if (!Stack.empty()) {
    mprinterr("Error: Mask stack is not empty.\n");
    while (!Stack.empty()) {
      delete[] Stack.top();
      Stack.pop();
    }
    delete[] pMask;
    return true;
  }

  if (intMask)
    maskIn.SetupMask( pMask, atoms_.size(), debug_ );
  else
    maskIn.SetupCharMask( pMask, atoms_.size(), debug_);
  delete[] pMask;
  return false;
}

// -----------------------------------------------------------------------------
// Topology::modifyStateByMask()
/**  The goal of this routine is to create a new AmberParm (newParm)
  *  based on the current AmberParm (this), deleting atoms that are
  *  not in the Selected array.
  */
Topology *Topology::modifyStateByMask(AtomMask &Mask) {
  std::vector<int> Map;

  //int newatom = 0;
  Map.reserve( Mask.Nselected() );
  for (AtomMask::const_iterator oldatom = Mask.begin(); oldatom != Mask.end(); oldatom++) 
    Map.push_back( *oldatom ); // Map[newatom] = oldatom

  return ModifyByMap(Map);
}

// Topology::ModifyByMap()
Topology *Topology::ModifyByMap(std::vector<int>& MapIn) {
  Topology *newParm = new Topology();

  newParm->parmName_ = parmName_;
  newParm->fileName_ = fileName_;
  newParm->radius_set_ = radius_set_;

  // Reverse Atom map
  // TODO: Use std::map instead
  std::vector<int> atomMap( atoms_.size(),-1 );

  // Copy atoms from this parm that are in Mask to newParm.
  //int newatom = 0;
  int oldres = -1;
  int oldmol = -1;
  newParm->firstSolventMol_ = -1;
  // TODO: Check the map size
  for (int newatom = 0; newatom < (int)MapIn.size(); newatom++) {
    int oldatom = MapIn[ newatom ];
    if (oldatom < 0) continue;
  //for (AtomMask::const_iterator oldatom = Mask.begin(); oldatom != Mask.end(); oldatom++) 
  //{
    // Store map of oldatom to newatom
    atomMap[oldatom] = newatom;
    // Copy oldatom 
    Atom newparmAtom = atoms_[oldatom];
    // Save oldatom residue number
    int curres = newparmAtom.ResNum();
    // Check if this old atom is in a different residue than the last. If so,
    // set new residue information.
    if ( curres != oldres ) {
      newParm->residues_.push_back( Residue( residues_[curres].Name(), newatom ) );
      oldres = curres;
    }
    // Clear bond information from new atom
    newparmAtom.ClearBonds();
    // Set new atom num and residue num
    //newparmAtom.SetNum( newatom );
    newparmAtom.SetResNum( newParm->residues_.size() - 1 );
    // Place new atom in newParm
    newParm->atoms_.push_back( newparmAtom );
    // Check if this old atom is in a different molecule than the last. If so,
    // set molecule information.
    int curmol = atoms_[oldatom].Mol();
    if (curmol != oldmol) {
      // Check if this is the first solvent mol of new parm
      if (newParm->firstSolventMol_==-1 && molecules_[curmol].IsSolvent()) {
        newParm->firstSolventMol_ = (int)newParm->molecules_.size();
        // Minus 2 since final solute residue is previous one and residues
        // has already been incremented. 
        newParm->finalSoluteRes_ = (int)newParm->residues_.size() - 2;
      }
      newParm->StartNewMol();
      oldmol = curmol;
    }
    // Copy extra amber info
    if (!itree_.empty()) newParm->itree_.push_back( itree_[oldatom] );
    if (!join_.empty()) newParm->join_.push_back( join_[oldatom] );
    if (!irotat_.empty()) newParm->irotat_.push_back( irotat_[oldatom] );
    //++newatom;
  }

  // NOTE: Since in the bond/angle/dihedral atom arrays the parm indices have 
  //       survived intact we can just include direct copies of all the 
  //       parameter arrays for now. May want to cull unused params later.

  // Set up new bond information
  newParm->bonds_ = SetupSequentialArray(atomMap, 3, bonds_);
  newParm->bondsh_ = SetupSequentialArray(atomMap, 3, bondsh_);
  newParm->SetAtomBondInfo();
  newParm->bondrk_ = bondrk_;
  newParm->bondreq_ = bondreq_;
  // Set new molecule information based on new bonds
  newParm->DetermineMolecules();
  // Set new solvent information based on new molecules
  newParm->SetSolventInfo(); 
  // Set up new angle info
  newParm->angles_ = SetupSequentialArray(atomMap, 4, angles_);
  newParm->anglesh_ = SetupSequentialArray(atomMap, 4, anglesh_);
  newParm->angletk_ = angletk_;
  newParm->angleteq_ = angleteq_;
  // Set up new dihedral info
  newParm->dihedrals_ = SetupSequentialArray(atomMap, 5, dihedrals_);
  newParm->dihedralsh_ = SetupSequentialArray(atomMap, 5, dihedralsh_);
  newParm->dihedralpk_ = dihedralpk_;
  newParm->dihedralpn_ = dihedralpn_;
  newParm->dihedralphase_ = dihedralphase_;
  newParm->scee_ = scee_;
  newParm->scnb_ = scnb_;
  // Set up nonbond info
  // Since nbindex depends on the atom type index and those entries were 
  // not changed this is still valid. May want to cull unused parms later.
  newParm->ntypes_ = ntypes_;
  newParm->nbindex_ = nbindex_;
  newParm->lja_ = lja_;
  newParm->ljb_ = ljb_;
  // Hbond info
  newParm->asol_ = asol_;
  newParm->bsol_ = bsol_;
  newParm->hbcut_ = hbcut_;
  // TODO: itree, join, irotat
  newParm->solty_ = solty_;
  
  // Setup excluded atoms list - Necessary?
  newParm->DetermineExcludedAtoms();

  // Give stripped parm the same pindex as original
  newParm->pindex_ = pindex_;

  newParm->nframes_ = nframes_;

  // Copy box information
  newParm->box_ = box_;

  return newParm;
}

// Topology::SetupSequentialArray()
/** Given an array with format [I0J0...X0][I1J1...] where the entries up to
  * X are atom# * 3 and entry X is an index, and an atom map array
  * with format Map[oldAtom]=newAtom, create a new array that contains
  * only entries for which all atoms are present. Can be used for the
  * bond, angle, and dihedral index arrays.
  */
std::vector<int> Topology::SetupSequentialArray(std::vector<int> &atomMap, int Nsequence,
                                                std::vector<int> &oldArray)
{
  std::vector<int> newArray;
  int Nsequence1 = Nsequence - 1;
  std::vector<int> newatoms(Nsequence, 0);
  std::vector<int> newsign(Nsequence1, 1);
  // Go through old array. Use atomMap to determine what goes into newArray.
  for (std::vector<int>::iterator oldi = oldArray.begin(); oldi != oldArray.end();
                                                           oldi += Nsequence)
  {
    // Using atomMap, check that atoms 0 to Nsequence exist in newParm. If
    // any of the atoms do not exist, bail.
    int newatm = -1;
    bool reverseOrder = false;
    std::vector<int>::iterator endidx = oldi + Nsequence1;
    unsigned int sequencei = 0;
    for (std::vector<int>::iterator oidx = oldi; oidx != endidx; oidx++) {
      int arrayIdx = *oidx;
      // For dihedrals the atom # can be negative. Convert to positive
      // for use in the atom map.
      if (arrayIdx < 0) {
        newatm = atomMap[ -arrayIdx / 3 ];
        newsign[sequencei] = -1;
        // For improper/multi-term dihedrals atom index 0 cannot be the 3rd
        // or 4th position since there is no such thing as -0.
        if (newatm == 0 && (sequencei == 2 || sequencei == 3))
          reverseOrder = true;
      } else {
        newatm = atomMap[ arrayIdx / 3 ];
        newsign[sequencei] = 1;
      }
      // New atom # of -1 means atom was removed - exit loop now.
      if (newatm == -1) break;
      // Atom exists - store.
      newatoms[sequencei++] = newatm;
    }
    // If newatm is -1 here that means it didnt exist in newParm for this
    // sequence. Skip the entire sequence.
    if (newatm==-1) continue;
    // Store the final number of the sequence, which is an index
    newatoms[Nsequence1] = *endidx;
    // Place the atoms in newatoms in newArray
    if (!reverseOrder) {
      for (int sequencei = 0; sequencei < Nsequence1; sequencei++)
        newArray.push_back( newatoms[sequencei] * 3 * newsign[sequencei] );
    } else {
      int ridx = Nsequence1 - 1;
      for (int sequencei = 0; sequencei < Nsequence1; sequencei++)
        newArray.push_back( newatoms[ridx--] * 3 * newsign[sequencei] );
    }
    // Place the index in newatoms in newArray
    newArray.push_back( newatoms[Nsequence1] );
  }
  return newArray;
}

