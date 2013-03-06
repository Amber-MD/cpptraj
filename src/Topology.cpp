#include <stack> // For ParseMask
#include <algorithm> // sort
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
#include "DistRoutines.h"
// DEBUG
//#include <cmath> // sqrt

// CONSTRUCTOR
Topology::Topology() :
  offset_(0.20),
  debug_(0),
  NsolventMolecules_(0),
  finalSoluteRes_(-1),
  pindex_(0),
  nframes_(0),
  ntypes_(0)
{ }

// Topology::SetOffset()
void Topology::SetOffset( double oIn ) {
  if (oIn > 0.0) offset_ = oIn;
}

// Topology::SetDebug()
void Topology::SetDebug(int debugIn) {
  debug_ = debugIn;
}

// Topology::SetParmName()
void Topology::SetParmName(std::string const& title, std::string const& filename) {
  parmName_ = title;
  fileName_ = filename;
}

// Topology::SetGBradiiSet()
void Topology::SetGBradiiSet(std::string const& gbset) {
  radius_set_ = gbset;
}  

// Topology::SetPindex()
void Topology::SetPindex(int pindexIn) {
  pindex_ = pindexIn;
}

// Topology::SetReferenceCoords()
void Topology::SetReferenceCoords( Frame* frameptr ) {
  if (frameptr == 0) return;
  refCoords_ = *frameptr;
}

// Topology::IncreaseFrames()
void Topology::IncreaseFrames(int frames) {
  nframes_ += frames;
}

// -----------------------------------------------------------------------------
// Topology::FinalSoluteRes()
/** Return 1 past the last solute residue. */
int Topology::FinalSoluteRes() const {
  return finalSoluteRes_ + 1;
}

// Topology::c_str()
/** Return a printf-compatible char* of the parm filename, or the parm
  * name (title) if the parm filename is empty.
  */
const char *Topology::c_str() const {
  if (!fileName_.empty()) 
    return fileName_.c_str();
  return parmName_.c_str();
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
  return GetBondLength( atoms_[atom1].Element(), atoms_[atom2].Element() );
}
// -----------------------------------------------------------------------------

// Topology::TruncResAtomName()
/** Given an atom number, return a string containing the corresponding 
  * residue name and number (starting from 1) along with the atom name 
  * with format: 
  * "<resname><resnum>@<atomname>", e.g. "ARG_11@CA".
  * Truncate the residue and atom names so there are no blanks.
  */
std::string Topology::TruncResAtomName(int atom) {
  std::string res_name;
  if (atom < 0 || atom >= (int)atoms_.size()) return res_name;
  // Atom name with no trailing spaces.
  std::string atom_name = atoms_[atom].Name().Truncated();
  int res = atoms_[atom].ResNum();
  // Residue name with no trailing spaces.
  // NOTE: ensure a residue size of 4?
  res_name = residues_[res].Name().Truncated();
  ++res; // want output as res+1
  res_name += "_";
  res_name += integerToString(res);
  res_name += "@";
  res_name += atom_name;
  return res_name;
}

// Topology::TruncResNameNum()
/** Given a residue number (starting from 0), return a string containing 
  * residue name and number (starting from 1) with format: 
  * "<resname>:<resnum>", e.g. "ARG:11".
  * Truncate residue name so there are no blanks.
  */
// FIXME: Add residue bounds check.
std::string Topology::TruncResNameNum(int res) {
  // Residue name with no trailing spaces.
  return residues_[res].Name().Truncated() + ":" + integerToString( res+1 );
}

// Topology::FindAtomInResidue()
/** Find the atom # of the specified atom name in the given residue.
  * \param res Residue number to search.
  * \param atname Atom name to find.
  * \return the atom number of the specified atom if found in the given residue.
  * \return -1 if atom not found in given residue.
  */
int Topology::FindAtomInResidue(int res, NameType const& atname) const {
  if (res < 0 || res >= (int)residues_.size()) return -1;
  for (int at = residues_[res].FirstAtom(); at < residues_[res].LastAtom(); ++at)
    if ( atoms_[at].Name() == atname )
      return at;
  return -1;
}

// Topology::FindResidueMaxNatom()
/** Return the # atoms in the largest residue. */
int Topology::FindResidueMaxNatom() const {
  if (residues_.size() <= 1)
    return (int)atoms_.size();
  int largest_natom = 0;
  for (std::vector<Residue>::const_iterator res = residues_.begin();
                                            res != residues_.end(); res++)
  {
    int diff = (*res).NumAtoms();
    if (diff > largest_natom) largest_natom = diff;
  }
  return largest_natom;
}

// Topology::SoluteAtoms()
// TODO: do not rely on finalsoluteres since it makes assumptions about
//       system layout
int Topology::SoluteAtoms() {
  if (NsolventMolecules_ == 0)
    return (int)atoms_.size();
  return ( residues_[finalSoluteRes_].LastAtom() );
}

// -----------------------------------------------------------------------------
// Topology::Summary()
void Topology::Summary() {
  mprintf("\t\tTopology %s contains %zu atoms.\n", c_str(), atoms_.size());
  mprintf("\t\t                  %zu residues.\n", residues_.size());
  mprintf("\t\t                  %zu bonds.\n", (bonds_.size()+bondsh_.size()) / 3 );
  mprintf("\t\t                  %zu molecules.\n", molecules_.size());
  mprintf("\t\t                  Box: %s\n",box_.TypeName());
  if (NsolventMolecules_>0) {
    mprintf("\t\t                  %i solvent molecules.\n", NsolventMolecules_);
    if (finalSoluteRes_>-1)
      mprintf("\t\t                  Final solute residue is %i\n", finalSoluteRes_+1);
  }
  /*if (!bondsh_.empty() || !bonds_.empty())
    mprintf("  %zu bonds to hydrogen, %zu other bonds.\n",bondsh_.size()/3,
            bonds_.size()/3);*/
}

// Topology::ParmInfo()
void Topology::ParmInfo() {
  mprintf(" %s, %zu atoms, %zu res, box: %s, %zu mol", c_str(),
          atoms_.size(), residues_.size(), box_.TypeName(), molecules_.size());
  if (NsolventMolecules_>0)
    mprintf(", %i solvent", NsolventMolecules_);
}

// Topology::PrintAtomInfo()
void Topology::PrintAtomInfo(std::string const& maskString) {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, true); // integer mask
  if ( mask.None() )
    mprintf("\tSelection is empty.\n");
  else {
    int width = DigitWidth(atoms_.size());
    if (width < 5) width = 5;
    mprintf("%-*s %4s %*s %4s %*s %4s %8s %8s\n", 
            width, "#Atom", "Name", 
            width, "#Res",  "Name",
            width, "#Mol",  "Type", "Charge", "Mass");
    for (AtomMask::const_iterator atnum = mask.begin(); atnum != mask.end(); atnum++) {
      const Atom& atom = atoms_[*atnum];
      int resnum = atom.ResNum();
      mprintf("%*i %4s %*i %4s %*i %4s %8.4f %8.4f\n", 
              width, *atnum+1, atom.c_str(), 
              width, resnum+1, residues_[resnum].c_str(),
              width, atom.Mol()+1, *(atom.Type()), atom.Charge(), atom.Mass());
    }
  }
}

// Topology::PrintBonds()
/** \param maskIn AtomMask which should have already been set up as a char mask
  */
void Topology::PrintBonds(std::vector<int>& barray, AtomMask const& maskIn) {
  for (std::vector<int>::iterator batom = barray.begin();
                                  batom != barray.end(); batom++)
  {
    int atom1 = ((*batom) / 3);
    if (!maskIn.AtomInCharMask( atom1 )) {
      batom += 2;
      continue;
    }
    ++batom;
    int atom2 = ((*batom) / 3);
    if (!maskIn.AtomInCharMask( atom2 )) {
      ++batom;
      continue;
    }
    mprintf("\tAtom %i:%s to %i:%s", atom1+1, atoms_[atom1].c_str(),
                                     atom2+1, atoms_[atom2].c_str());
    ++batom;
    if (*batom==-1) {
      double req = GetBondLength(atoms_[atom1].Element(),atoms_[atom2].Element());
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
void Topology::PrintBondInfo(std::string const& maskString) {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, false); // Char mask
  if (!bondsh_.empty()) {
    mprintf("%zu BONDS TO HYDROGEN:\n",bondsh_.size()/3);
    PrintBonds( bondsh_, mask );
  }
  if (!bonds_.empty()) {
    mprintf("%zu BONDS TO NON-HYDROGEN:\n",bonds_.size()/3);
    PrintBonds( bonds_, mask );
  }
}

// Topology::PrintMoleculeInfo()
void Topology::PrintMoleculeInfo(std::string const& maskString) {
  if (molecules_.empty())
    mprintf("\t[%s] No molecule info.\n",c_str());
  else {
    AtomMask mask( maskString );
    ParseMask(refCoords_, mask, false); // Char mask
    mprintf("MOLECULES:\n");
    unsigned int mnum = 1;
    for (std::vector<Molecule>::iterator mol = molecules_.begin(); 
                                         mol != molecules_.end(); mol++)
    {
      if ( mask.AtomsInCharMask( (*mol).BeginAtom(), (*mol).EndAtom() ) ) {
        int firstres = atoms_[ (*mol).BeginAtom() ].ResNum();
        mprintf("\tMolecule %u, %i atoms, first residue %i:%s",mnum,(*mol).NumAtoms(),
                firstres+1, residues_[firstres].c_str());
        if ( (*mol).IsSolvent() ) mprintf(" SOLVENT");
        mprintf("\n");
      }
      ++mnum;
    }
  }
}

// Topology::PrintResidueInfo()
void Topology::PrintResidueInfo(std::string const& maskString) {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, false); // Char mask
  mprintf("RESIDUES:\n");
  unsigned int rnum = 1;
  for (std::vector<Residue>::iterator res = residues_.begin();
                                      res != residues_.end(); res++)
  {
    if ( mask.AtomsInCharMask( (*res).FirstAtom(), (*res).LastAtom() ) ) {
      mprintf("\tResidue %u %s first atom %i last atom %i\n",
              rnum, (*res).c_str(), (*res).FirstAtom()+1, (*res).LastAtom());
    }
    ++rnum;
  }
}

void Topology::PrintChargeInfo(std::string const& maskString) {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, true); // Int mask
  double sumq = 0.0;
  for (AtomMask::const_iterator aidx = mask.begin(); aidx != mask.end(); ++aidx)
    sumq += atoms_[*aidx].Charge();
  mprintf("\tSum of charges in mask");
  mask.BriefMaskInfo();
  mprintf(" is %f\n", sumq);
}

// -----------------------------------------------------------------------------
// Topology::AddTopAtom()
void Topology::AddTopAtom(Atom atomIn, NameType const& resname, int current_res, int& last_res, 
                       const double* XYZin) 
{
  // Check if this is a new residue
  if (residues_.empty() || current_res != last_res) {
    // Last atom of old residue is == current # atoms.
    if (!residues_.empty())
      residues_.back().SetLastAtom( atoms_.size() );
    // First atom of new residue is == current # atoms.
    residues_.push_back( Residue(resname, atoms_.size()) );
    last_res = current_res;
  }
  // Set this atoms residue number 
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
    molecules_.push_back( Molecule(0, atoms_.size()) );
  } else {
    // The first atom of this molecule will be end atom of last molecule.
    int molFirstAtom = molecules_.back().EndAtom();
    // Only add a new molecule if #atoms > first atom of the molecule.
    if ((int)atoms_.size() > molFirstAtom) 
      molecules_.push_back( Molecule( molFirstAtom, atoms_.size()) );
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
  for (size_t res = 0; res < nres - 1; res++) 
    residues_.push_back( Residue( resnames[res], resnums[res], resnums[res+1] ) );
  residues_.push_back( Residue( resnames[nres-1], resnums[nres-1], atoms_.size() ) );

  return 0;
}

// Topology::SetBondInfo()
int Topology::SetBondInfo(std::vector<int> &bonds, std::vector<int> &bondsh,
                          std::vector<double> &bond_rk, std::vector<double> &bond_req) 
{
  if (bonds.empty() && bondsh.empty())
    mprinterr("Warning: Topology: Input bonds and bondsh are empty.\n");

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
int Topology::CommonSetup(bool bondsearch) {
  // Set residue last atom (PDB/Mol2/PSF) 
  residues_.back().SetLastAtom( atoms_.size() );
  // Set up bond information if specified and necessary
  if (bondsearch) {
    if (bonds_.empty() && bondsh_.empty() && !refCoords_.empty()) {
      GetBondsFromAtomCoords();
      if (DetermineMolecules()) return 1;
    }
  }
  
  if (molecules_.empty()) 
    if (DetermineMolecules()) return 1;

  // Set up solvent information
  if (SetSolventInfo()) return 1;

  // Determine excluded atoms
  DetermineExcludedAtoms();

  return 0;
}

// -----------------------------------------------------------------------------
// WarnBondLengthDefault()
void Topology::WarnBondLengthDefault(Atom::AtomicElementType atom1,
                                     Atom::AtomicElementType atom2, double cut) {
  mprintf("Warning: GetBondLength: Bond length not found for %s - %s\n",
          Atom::AtomicElementName[atom1], Atom::AtomicElementName[atom2]);
  mprintf("                        Using default length of %f\n", cut);
}

// Topology::GetBondLength() 
/** Return optimal covalent bond distance based on the element types of atom1 
  * and atom2. When multiple hybridizations are possible the longest possible 
  * bond length is used.
  * Unless otherwise noted values taken from:
  * - Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 
  *       2nd ed., Butterworths, London, 1958; 
  * - B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, 
  *       No. 31, Washington, DC, 1970; S.W. Benson, J. Chem. Educ., 42, 502 (1965).
  * Can be found on the web at:
  * - http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
  */
// NOTE: Store cut^2 instead??
double Topology::GetBondLength(Atom::AtomicElementType atom1, Atom::AtomicElementType atom2) 
{
  // Default cutoff
  double cut = 1.60;
  if (atom1==atom2) {
    // Self
    switch (atom1) {
      case Atom::HYDROGEN  : cut=0.74; break;
      case Atom::CARBON    : cut=1.54; break;
      case Atom::NITROGEN  : cut=1.45; break;
      case Atom::OXYGEN    : cut=1.48; break;
      case Atom::PHOSPHORUS: cut=2.21; break;
      case Atom::SULFUR    : cut=2.05; break; // S-S gas-phase value; S=S is 1.49
      default: WarnBondLengthDefault(atom1,atom2,cut);
    }
  } else {
    Atom::AtomicElementType e1, e2;
    if (atom1 < atom2) {
      e1 = atom1;
      e2 = atom2;
    } else {
      e1 = atom2;
      e2 = atom1;
    }
    switch (e1) {
      case Atom::HYDROGEN: // Bonds to H
        switch (e2) {
          case Atom::CARBON    : cut=1.09; break;
          case Atom::NITROGEN  : cut=1.01; break;
          case Atom::OXYGEN    : cut=0.96; break;
          case Atom::PHOSPHORUS: cut=1.44; break;
          case Atom::SULFUR    : cut=1.34; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case Atom::CARBON: // Bonds to C
        switch (e2) {
          case Atom::NITROGEN  : cut=1.47; break;
          case Atom::OXYGEN    : cut=1.43; break;
          case Atom::FLUORINE  : cut=1.35; break;
          case Atom::PHOSPHORUS: cut=1.84; break;
          case Atom::SULFUR    : cut=1.82; break;
          case Atom::CHLORINE  : cut=1.77; break;
          case Atom::BROMINE   : cut=1.94; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case Atom::NITROGEN: // Bonds to N
        switch (e2) {
          case Atom::OXYGEN    : cut=1.40; break;
          case Atom::FLUORINE  : cut=1.36; break;
          case Atom::PHOSPHORUS: cut=1.71; // Avg over all nX-pX from gaff.dat
          case Atom::SULFUR    : cut=1.68; break; // Postma & Vos, Acta Cryst. (1973) B29, 915
          case Atom::CHLORINE  : cut=1.75; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case Atom::OXYGEN: // Bonds to O
        switch (e2) {
          case Atom::FLUORINE  : cut=1.42; break;
          case Atom::PHOSPHORUS: cut=1.63; break;
          case Atom::SULFUR    : cut=1.48; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case Atom::FLUORINE: // Bonds to F
        switch (e2) {
          case Atom::PHOSPHORUS: cut=1.54; break;
          case Atom::SULFUR    : cut=1.56; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;  
      case Atom::PHOSPHORUS: // Bonds to P
        switch (e2) {
          case Atom::SULFUR  : cut=1.86; break;
          case Atom::CHLORINE: cut=2.03; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case Atom::SULFUR: // Bonds to S
        switch (e2) {
          case Atom::CHLORINE: cut=2.07; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      default: WarnBondLengthDefault(e1,e2,cut);
    } // END switch(e1)
  }
  //mprintf("\t\tCUTOFF: [%s] -- [%s] = %lf\n",AtomicElementName[atom1],
  //        AtomicElementName[atom2],cut);
  return cut;
}

bool Topology::NameIsSolvent(NameType const& resname) {
  if (resname == "WAT " ||
      resname == "HOH " ||
      resname == "TIP3"
     )
    return true;
  return false;
}

// Topology::GetBondsFromAtomCoords()
void Topology::GetBondsFromAtomCoords() {
  mprintf("\t%s: determining bond info from distances.\n",c_str());
  // ----- STEP 1: Determine bonds within residues
  for (std::vector<Residue>::iterator res = residues_.begin(); 
                                      res != residues_.end(); ++res) 
  {
    // Get residue start atom.
    int startatom = (*res).FirstAtom();
    // Get residue end atom.
    int stopatom = (*res).LastAtom();
    // DEBUG
    //mprintf("\tRes %i Start atom %zu coords: ",resnum+1, startatom+1);
    //atoms_[startatom].PrintXYZ();
    //mprintf("\n");
    // Check for bonds between each atom in the residue.
    for (int atom1 = startatom; atom1 < stopatom - 1; ++atom1) {
      // If this is a hydrogen and it already has a bond, move on.
      if (atoms_[atom1].Element()==Atom::HYDROGEN &&
          atoms_[atom1].Nbonds() > 0 )
        continue;
      for (int atom2 = atom1 + 1; atom2 < stopatom; ++atom2) {
        double D2 = DIST2_NoImage(refCoords_.XYZ(atom1), refCoords_.XYZ(atom2) );
        double cutoff2 = GetBondLength(atoms_[atom1].Element(), atoms_[atom2].Element()) + offset_;
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
  }

  // ----- STEP 2: Determine bonds between adjacent residues
  std::vector<Molecule>::iterator nextmol = molecules_.begin();
  if (!molecules_.empty())
    ++nextmol;
  for (std::vector<Residue>::iterator res = residues_.begin() + 1;
                                      res != residues_.end(); ++res)
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
    // If this residue is recognized as solvent, no need to check previous or
    // next residue
    if ( NameIsSolvent((*res).Name()) ) {
      ++res;
      if (res == residues_.end()) break;
      continue;
    }
    // Get previous residue
    std::vector<Residue>::iterator previous_res = res - 1;
    // If previous residue is recognized as solvent, no need to check previous.
    if ( NameIsSolvent((*previous_res).Name()) ) continue;
    // Get previous residue start atom
    int startatom = (*previous_res).FirstAtom();
    // Previous residue stop atom, this residue start atom
    int midatom = (*res).FirstAtom();
    // This residue stop atom
    int stopatom = (*res).LastAtom();
    //mprintf("\tBonds between residues %s and %s\n",(*previous_res).c_str(),(*res).c_str());
    // Check for bonds between adjacent residues
    for (int atom1 = startatom; atom1 < midatom; atom1++) {
      if (atoms_[atom1].Element()==Atom::HYDROGEN) continue;
      for (int atom2 = midatom; atom2 < stopatom; atom2++) {
        if (atoms_[atom2].Element()==Atom::HYDROGEN) continue;
        double D2 = DIST2_NoImage(refCoords_.XYZ(atom1), refCoords_.XYZ(atom2) );
        double cutoff2 = GetBondLength(atoms_[atom1].Element(), atoms_[atom2].Element()) + offset_;
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
int Topology::DetermineMolecules() {
  std::vector<Atom>::iterator atom;

  if (debug_ > 0) mprintf("\t%s: determining molecule info from bonds.\n",c_str());
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
  if (mol == 0) return 0;
  std::vector<Molecule>::iterator molecule = molecules_.begin();
  (*molecule).SetFirst(0);
  atom = atoms_.begin(); 
  int lastMol = (*atom).Mol();
  int atomNum = 0;
  for (; atom != atoms_.end(); atom++)
  {
    if ( (*atom).Mol() > lastMol ) {
      // Set last atom of molecule
      (*molecule).SetLast( atomNum );
      // Set first atom of next molecule
      ++molecule;
      (*molecule).SetFirst( atomNum );
      lastMol = (*atom).Mol();
    } else if ( (*atom).Mol()  < lastMol) {
      mprinterr("Error: Atom %u was assigned a lower molecule # than previous atom.\n",
                atom - atoms_.begin() + 1);
      mprinterr("Error: This can happen if bond information is incorrect or missing.\n");
      mprinterr("Error: Detected # of molecules is %i. If this is incorrect and your\n",
                mol);
      mprinterr("Error: topology does not have bond information (e.g. PDB file), try\n");
      mprinterr("Error: increasing the bond search cutoff offset (currently %f).\n",offset_);
      mprinterr("Error: e.g. 'parm %s bondsearch <new offset>'\n", fileName_.c_str());
      molecules_.clear();
      return 1;
    }
    ++atomNum;
  }
  (*molecule).SetLast( atoms_.size() );
  return 0;
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
// Topology::SetSolvent()
/** Set solvent information from atom mask. */
int Topology::SetSolvent(std::string const& maskexpr) {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolvent [%s]: No molecule information.\n", c_str());
    return 1;
  }
  // Setup mask
  AtomMask mask( maskexpr );
  SetupCharMask( mask );
  if (mask.None()) {
    mprinterr("Error: SetSolvent [%s]: Mask %s selects no atoms.\n", c_str(), maskexpr.c_str());
    return 1;
  }
  // Loop over all molecules
  NsolventMolecules_ = 0;
  finalSoluteRes_ = -1;
  int numSolvAtoms = 0;
  int firstSolventMol = -1;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); ++mol)
  {
    // Reset old solvent information.
    (*mol).SetNoSolvent();
    // If any atoms in this molecule are selected by mask, make entire
    // molecule solvent.
    for (int atom = (*mol).BeginAtom(); atom < (*mol).EndAtom(); ++atom) {
      if ( mask.AtomInCharMask( atom ) ) {
        (*mol).SetSolvent();
        if (firstSolventMol == -1) {
          // This is first solvent mol. Final solute res is the one before this.
          int firstRes = atoms_[ (*mol).BeginAtom() ].ResNum();
          finalSoluteRes_ = firstRes - 1;
          firstSolventMol = (int)(mol - molecules_.begin());
        }
        ++NsolventMolecules_;
        numSolvAtoms += (*mol).NumAtoms();
        break;
      }
    }
  }

  if (firstSolventMol == -1 && finalSoluteRes_ == -1)
    finalSoluteRes_ = (int)residues_.size() - 1;
  mprintf("\tSolvent Mask [%s]: %i solvent molecules, %i solvent atoms\n",
          maskexpr.c_str(), NsolventMolecules_, numSolvAtoms);
  return 0;
}

// Topology::SetSolventInfo()
/** Determine which molecules are solvent based on residue name. 
  * Also set finalSoluteRes.
  */
int Topology::SetSolventInfo() {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolventInfo: No molecule information.\n");
    return 1;
  }
  // Loop over each molecule. Check if first residue of molecule is solvent.
  NsolventMolecules_ = 0;
  finalSoluteRes_ = -1;
  int numSolvAtoms = 0;
  int firstSolventMol = -1;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); mol++)
  {
    int firstRes = atoms_[ (*mol).BeginAtom() ].ResNum();
    if ( NameIsSolvent(residues_[firstRes].Name()) ) {
      (*mol).SetSolvent();
      ++NsolventMolecules_;
      numSolvAtoms += (*mol).NumAtoms();
      if (firstSolventMol==-1) {
        // This is first solvent mol. Final solute res is the one before this.
        finalSoluteRes_ = firstRes - 1;
        firstSolventMol = (int)(mol - molecules_.begin());
      }
    }
  }

  if (firstSolventMol == -1 && finalSoluteRes_ == -1)
    finalSoluteRes_ = (int)residues_.size() - 1;
  if (debug_>0) {
    if (NsolventMolecules_ == 0) 
      mprintf("    No solvent.\n");
    else
      mprintf("    %i solvent molecules, %i solvent atoms\n",NsolventMolecules_,numSolvAtoms);
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask) const { 
  return ParseMask(refCoords_, mask, true);
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask) const {
  return ParseMask(refCoords_, mask, false);
}

// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask, Frame const& frame) const {
  return ParseMask( frame, mask, true );
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask, Frame const& frame) const {
  return ParseMask( frame, mask, false );
}

// Topology::Mask_SelectDistance()
void Topology::Mask_SelectDistance( Frame const& REF, char *mask, bool within, 
                                    bool byAtom, double distance ) const 
{
  int endatom, resi;
  bool selectresidue;
  int atomi, idx, atomj;
  double d2;
  const double* i_crd;

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
    int n_of_atoms = (int)atoms_.size();
#ifdef _OPENMP
#pragma omp parallel private(atomi, idx, atomj, d2, i_crd)
{
#pragma omp for
#endif
    for (atomi = 0; atomi < n_of_atoms; atomi++) {
      // No need to calculate if atomi already selected
      if (mask[atomi] == 'T') continue;
      // Loop over initially selected atoms
      i_crd = REF.XYZ( atomi );
      for (idx = 0; idx < (int)selected.size(); idx++) {
        atomj = selected[idx];
        d2 = DIST2_NoImage(i_crd, REF.XYZ(atomj));
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
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  } else { // Select by residue
    int n_of_res = (int)residues_.size();
#ifdef _OPENMP
#pragma omp parallel private(atomi, idx, atomj, d2, resi, selectresidue, endatom, i_crd)
{
#pragma omp for
#endif
    for (resi = 0; resi < n_of_res; resi++) {
      selectresidue = false;
      // Determine end atom for this residue
      endatom = residues_[resi].LastAtom();
      // Loop over mask atoms
      for (idx = 0; idx < (int)selected.size(); idx++) {
        atomj = selected[idx];
        i_crd = REF.XYZ( atomj );
        // Loop over residue atoms
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++) {
          d2 = DIST2_NoImage(REF.XYZ(atomi), i_crd);
          if (within) {
            if (d2 < distance) selectresidue = true;
          } else {
            if (d2 > distance) selectresidue = true;
          }
          if (selectresidue) break; 
        }
        if (selectresidue) break;
      }
      if (selectresidue) {
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++)
          mask[atomi] = 'T';
        continue;
      }
    }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  }
}

// Topology::Mask_AND()
void Topology::Mask_AND(char *mask1, char *mask2) const {
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
void Topology::Mask_OR(char *mask1, char *mask2) const {
  //mprintf("\t\t\tPerforming OR on masks.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T' || mask2[i]=='T')
      mask1[i] = 'T';
    else
      mask1[i] = 'F';
  }
}

// Topology::Mask_NEG()
void Topology::Mask_NEG(char *mask1) const {
  //mprintf("\t\t\tNegating mask.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T')
      mask1[i] = 'F';
    else
      mask1[i] = 'T';
  }
}

// Topology::MaskSelectResidues()
void Topology::MaskSelectResidues(NameType const& name, char *mask) const {
  //mprintf("\t\t\tSelecting residues named [%s]\n",*name);
  for (std::vector<Residue>::const_iterator res = residues_.begin();
                                            res != residues_.end(); res++)
  {
    if ( (*res).Name().Match( name ) ) {
      std::fill(mask + (*res).FirstAtom(), mask + (*res).LastAtom(), 'T');
    }
  }
}

// Topology::MaskSelectResidues()
// Mask args expected to start from 1
void Topology::MaskSelectResidues(int res1, int res2, char *mask) const {
  int endatom;
  int nres = (int) residues_.size();
  //mprintf("\t\t\tSelecting residues %i to %i\n",res1,res2);
  // Check start atom. res1 and res2 are checked by MaskToken
  if (res1 > nres) {
    if (debug_>0)
      mprintf("Warning: Select residues: res 1 out of range (%i)\n",res1);
    return;
  }
  // If last res > nres, make it nres
  if ( res2 >= nres )
    endatom = (int)atoms_.size();
  else
    endatom = residues_[res2-1].LastAtom();
  // Select atoms
  std::fill(mask + residues_[res1-1].FirstAtom(), mask + endatom, 'T');
}

// Topology::MaskSelectElements()
void Topology::MaskSelectElements( NameType const& element, char* mask ) const {
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
  {
    NameType atom_element( (*atom).ElementName() );
    if ( atom_element.Match( element ) )
      mask[m] = 'T';
    ++m;
  } 
}

// Topology::MaskSelectTypes()
void Topology::MaskSelectTypes( NameType const& type, char* mask ) const {
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
  {
    if ( (*atom).Type().Match( type ) )
      mask[m] = 'T';
    ++m;
  } 
}

// Topology::MaskSelectAtoms()
void Topology::MaskSelectAtoms( NameType const& name, char *mask) const {
  //mprintf("\t\t\tSelecting atoms named [%s]\n",*name);
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
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
void Topology::MaskSelectAtoms(int atom1, int atom2, char *mask) const {
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
bool Topology::ParseMask(Frame const& REF, AtomMask &maskIn, bool intMask) const {
  std::stack<char*> Stack;
  char *pMask = 0; 
  char *pMask2 = 0;

  for (AtomMask::token_iterator token = maskIn.begintoken();
                                token != maskIn.endtoken(); token++)
  {
    if (pMask==0) {
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
      case MaskToken::AtomType :
        MaskSelectTypes( (*token).Name(), pMask );
        break;
      case MaskToken::AtomElement :
        MaskSelectElements( (*token).Name(), pMask );
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
        mprinterr("Error: Invalid mask token (Mask [%s], type [%s]).\n",
                  maskIn.MaskString(), (*token).TypeName() );
    }
    // Check if this mask should now go on the stack
    if ( (*token).OnStack() ) {
      //mprintf("Pushing Mask on stack, last Token [%s]\n",(*token).TypeName());
      Stack.push( pMask );
      pMask = 0;
    }
  }
  // If pMask is not null it is probably a blank leftover
  if (pMask!=0) delete[] pMask;

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
    maskIn.SetupIntMask( pMask, atoms_.size(), debug_ );
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
  int firstSolventMol = -1;
  // TODO: Check the map size
  for (int newatom = 0; newatom < (int)MapIn.size(); newatom++) {
    int oldatom = MapIn[ newatom ];
    if (oldatom < 0) continue;
    // Store map of oldatom to newatom
    atomMap[oldatom] = newatom;
    // Copy oldatom 
    Atom newparmAtom = atoms_[oldatom];
    // Save oldatom residue number
    int curres = newparmAtom.ResNum();
    // Check if this old atom is in a different residue than the last. If so,
    // set new residue information.
    if ( curres != oldres ) {
      if (!newParm->residues_.empty())
        newParm->residues_.back().SetLastAtom( newatom );
      newParm->residues_.push_back( Residue(residues_[curres].Name(), newatom) );
      oldres = curres;
    }
    // Clear bond information from new atom
    newparmAtom.ClearBonds();
    // Set new atom num and residue num
    newparmAtom.SetResNum( newParm->residues_.size() - 1 );
    // Place new atom in newParm
    newParm->atoms_.push_back( newparmAtom );
    // Check if this old atom is in a different molecule than the last. If so,
    // set molecule information.
    int curmol = atoms_[oldatom].Mol();
    if (curmol != oldmol) {
      // Check if this is the first solvent mol of new parm
      if (firstSolventMol==-1 && molecules_[curmol].IsSolvent()) {
        firstSolventMol = (int)newParm->molecules_.size();
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
  }
  // Set last residue last atom
  newParm->residues_.back().SetLastAtom( newParm->atoms_.size() );

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
  if (newParm->DetermineMolecules()) {
    delete newParm;
    return 0;
  }
  // Set new solvent information based on new molecules
  if (newParm->SetSolventInfo()) {
    delete newParm;
    return 0;
  } 
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
  // NOTE: SOLTY is currently unused 
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

