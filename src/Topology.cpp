#include <cmath> // pow
#include <algorithm> // find
#include <map> // For atom types in append
#include <stack> // For large system molecule search
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
#include "Constants.h" // RADDEG, SMALL
#include "AtomType.h"

const NonbondType Topology::LJ_EMPTY = NonbondType();

// CONSTRUCTOR
Topology::Topology() :
  debug_(0),
  ipol_(0),
  NsolventMolecules_(0),
  pindex_(0),
  n_extra_pts_(0),
  n_atom_types_(0)
{ }

// Topology::SetParmName()
void Topology::SetParmName(std::string const& title, FileName const& filename) {
  parmName_ = title;
  fileName_ = filename;
}

// Topology::c_str()
/** Return a printf-compatible char* of the parm filename, or the parm
  * name (title) if the parm filename is empty.
  */
const char *Topology::c_str() const {
  if (!fileName_.empty())
    return fileName_.base();
  return parmName_.c_str();
}

/** Used to set box info from currently associated trajectory. */
// FIXME: This routine is here for potential backwards compatibility issues
//        since the topology box information was previously modified by
//        trajectory box information, but may no longer be necessary or
//        desirable.
void Topology::SetBoxFromTraj(Box const& boxIn) {
  if (!boxIn.HasBox()) {
    // No incoming box.
    if ( parmBox_.HasBox()) {
      // No incoming box and parm has box - disable parm box.
      mprintf("Warning: Box information present in topology but not in trajectory.\n"
              "Warning: DISABLING BOX in topology '%s'!\n", c_str());
      parmBox_.SetNoBox();
    }
  } else {
    // Incoming box.
    if ( boxIn.BoxX() < Constants::SMALL || 
         boxIn.BoxY() < Constants::SMALL || 
         boxIn.BoxZ() < Constants::SMALL )
    {
      // Incoming box has no lengths - disable parm box.
      mprintf("Warning: Box information present in trajectory but lengths are zero.\n"
              "Warning: DISABLING BOX in topology '%s'!\n", c_str());
      parmBox_.SetNoBox();
    } else {
      // Incoming box is valid. Indicate if current box type differs from
      // incoming box type.
      if (parmBox_.Type() != boxIn.Type()) {
        mprintf("Warning: Trajectory box type is '%s' but topology box type is '%s'.\n"
                "Warning: Setting topology box information from trajectory.\n",
                boxIn.TypeName(), parmBox_.TypeName());
      }
      parmBox_ = boxIn;
    }
  }
}

// Topology::SetDistMaskRef()
void Topology::SetDistMaskRef( Frame const& frameIn ) {
  if (!frameIn.empty()) {
    if (frameIn.Natom() == Natom())
      refCoords_ = frameIn;
    else if (frameIn.Natom() > Natom()) {
      mprintf("Warning: Active reference has %i atoms, parm '%s' has only %i.\n"
              "Warning: Truncating reference coords for this parm (distance-based masks only).\n",
              frameIn.Natom(), c_str(), Natom());
      refCoords_.SetupFrame(Natom());
      std::copy(frameIn.xAddress(), frameIn.xAddress() + refCoords_.size(),
                refCoords_.xAddress());
    } else {
      mprintf("Warning: Active reference has only %i atoms, parm '%s' has %i.\n"
              "Warning: Parm will only have reference coordinates for the first %i atoms"
              " (distance-based masks only).\n",
              frameIn.Natom(), c_str(), Natom(), frameIn.Natom());
      refCoords_.SetupFrame(Natom());
      std::copy(frameIn.xAddress(), frameIn.xAddress() + frameIn.size(), refCoords_.xAddress());
      std::fill(refCoords_.xAddress() + frameIn.size(),
                refCoords_.xAddress() + refCoords_.size(), 0.0);
    }
  }
}

// -----------------------------------------------------------------------------
/** \return Range containing only solute residues. */
Range Topology::SoluteResidues() const {
  Range solute_res;
  atom_iterator atom = atoms_.begin();
  while (atom != atoms_.end()) {
    // If atom is in a solvent molecule skip molecule. Otherwise add res num
    // and skip to next residue.
    if (molecules_[atom->MolNum()].IsSolvent())
      atom += molecules_[atom->MolNum()].NumAtoms();
    else if (molecules_[atom->MolNum()].NumAtoms() == 1) // Assume ion.
      ++atom;
    else {
      solute_res.AddToRange( atom->ResNum() );
      if (debug_ > 0)
        mprintf("DEBUG:\t\tAdding solute residue %i\n", atom->ResNum()+1);
      atom += residues_[atom->ResNum()].NumAtoms();
    }
  }
  return solute_res;
}

// -----------------------------------------------------------------------------
// Topology::TruncResAtomName()
/** Given an atom number, return a string containing the corresponding 
  * residue name and number (starting from 1) along with the atom name 
  * with format: 
  * "<resname>_<resnum>@<atomname>", e.g. "ARG_11@CA".
  * Truncate the residue and atom names so there are no blanks.
  */
std::string Topology::TruncResAtomName(int atom) const {
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

// Topology::TruncResAtomNameNum()
/** Given an atom number, return a string containing the corresponding 
  * residue name and number (starting from 1) along with the atom name 
  * and number with format: 
  * "<resname>_<resnum>@<atomname>_<atomnum>", e.g. "ARG_11@CA_256".
  * Truncate the residue and atom names so there are no blanks.
  */
std::string Topology::TruncResAtomNameNum(int atom) const {
  return TruncResAtomName(atom) + "_" + integerToString(atom+1);
}

// Topology::AtomMaskName()
/** \return A string of format :r@a where r is atoms residue number and
  *         a is atoms name.
  */
std::string Topology::AtomMaskName(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string(""); 
  std::string maskName = ":";
  maskName += integerToString( atoms_[atom].ResNum() + 1 );
  maskName += "@";
  maskName += atoms_[atom].Name().Truncated();
  return maskName;
}

std::string Topology::TruncAtomNameNum(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string("");
  std::string atom_name = atoms_[atom].Name().Truncated();
  atom_name += "_";
  atom_name += integerToString(atom + 1);
  return atom_name;
}

// Topology::TruncResNameNum()
/** Given a residue number (starting from 0), return a string containing 
  * residue name and number (starting from 1) with format: 
  * "<resname>:<resnum>", e.g. "ARG:11".
  * Truncate residue name so there are no blanks.
  */
// FIXME: Add residue bounds check.
std::string Topology::TruncResNameNum(int res) const {
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

// -----------------------------------------------------------------------------
// Topology::Summary()
void Topology::Summary() const {
  mprintf("\tTopology %s contains %zu atoms.\n", c_str(), atoms_.size());
  if (!parmName_.empty())
    mprintf("\t\tTitle: %s\n", parmName_.c_str());
  if (!fileName_.empty())
    mprintf("\t\tOriginal filename: %s\n", fileName_.full());
  mprintf("\t\t%zu residues.\n", residues_.size());
  mprintf("\t\t%zu molecules.\n", molecules_.size());
  size_t s1 = bondsh_.size();
  size_t s2 = bonds_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu bonds (%zu to H, %zu other).\n", s1+s2, s1, s2);
  s1 = anglesh_.size();
  s2 = angles_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu angles (%zu with H, %zu other).\n", s1+s2, s1 ,s2);
  s1 = dihedralsh_.size();
  s2 = dihedrals_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu dihedrals (%zu with H, %zu other).\n", s1+s2, s1, s2);
  mprintf("\t\tBox: %s\n", parmBox_.TypeName());
  if (NsolventMolecules_>0) {
    mprintf("\t\t%i solvent molecules.\n", NsolventMolecules_);
  }
  if (!radius_set_.empty())
    mprintf("\t\tGB radii set: %s\n", radius_set_.c_str());
  if (chamber_.HasChamber()) {
    mprintf("\t\tCHAMBER: %zu Urey-Bradley terms, %zu Impropers\n",
            chamber_.UB().size(), chamber_.Impropers().size());
    if (chamber_.HasCmap())
      mprintf("\t\t         %zu CMAP grids, %zu CMAP terms.\n", 
              chamber_.CmapGrid().size(), chamber_.Cmap().size());
  }
  if (lesparm_.HasLES())
    mprintf("\t\tLES info: %i types, %i copies\n", lesparm_.Ntypes(), lesparm_.Ncopies());
  if (cap_.HasWaterCap())
    mprintf("\t\tCAP info: Last atom before cap = %s, Cut= %g, X= %g, Y= %g, Z= %g\n",
            AtomMaskName(cap_.NatCap()).c_str(), cap_.CutCap(), 
            cap_.xCap(), cap_.yCap(), cap_.zCap());
}

// Topology::Brief()
void Topology::Brief(const char* heading) const {
  if (heading != 0)
    mprintf("\t%s", heading);
  else
    mprintf(" %s,", c_str());
  mprintf(" %zu atoms, %zu res, box: %s, %zu mol", atoms_.size(), 
          residues_.size(), parmBox_.TypeName(), molecules_.size());
  if (NsolventMolecules_>0)
    mprintf(", %i solvent", NsolventMolecules_);
  if (heading != 0)
    mprintf("\n");
}

// -----------------------------------------------------------------------------
// Topology::AddTopAtom()
int Topology::AddTopAtom(Atom const& atomIn, Residue const& resIn)
{
  // If no residues or res num has changed, this is a new residue.
  // TODO check chain ID?
  if ( residues_.empty() || 
       residues_.back().OriginalResNum() != resIn.OriginalResNum() ||
       residues_.back().SegID() != resIn.SegID() ||
       residues_.back().Icode() != resIn.Icode() )
  {
    // Last atom of old residue is == current # atoms.
    if (!residues_.empty())
      residues_.back().SetLastAtom( atoms_.size() );
    // First atom of new residue is == current # atoms.
    residues_.push_back( resIn );
    residues_.back().SetFirstAtom( atoms_.size() );
  }
  atoms_.push_back(atomIn);
  // Set this atoms internal residue number 
  atoms_.back().SetResNum( residues_.size()-1 );
  // Set current residues final atom number
  residues_.back().SetLastAtom( atoms_.size() );
  return 0;
}

// Topology::StartNewMol()
void Topology::StartNewMol() {
  // No atoms, so no need to do anything.
  if (atoms_.empty()) return;
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
  if (residues_.empty()) {
    // No residues yet. Consider entire molecule to be the residue.
    mprintf("Warning: Starting a molecule before residue info present.\n"
            "Warning:   Creating residue named 'MOL'\n");
    residues_.push_back( Residue("MOL",0,atoms_.size(),1,' ',' ') );
  } 
  residues_.back().SetTerminal( true );
}

// Topology::CommonSetup()
int Topology::CommonSetup(bool molsearch) {
  // TODO: Make bond parm assignment / molecule search optional?
  // Assign default lengths if necessary (for e.g. CheckStructure)
  if (bondparm_.empty())
    AssignBondParameters();
  if (molsearch) {
    // Determine molecule info from bonds
    if (DetermineMolecules())
      mprinterr("Error: Could not determine molecule information for %s.\n", c_str());
  }
  // Check that molecules do not share residue numbers. Only when bond searching.
  // FIXME always check? 
  if (!molecules_.empty() && molecules_.size() > 1) {
    bool mols_share_residues = (molecules_.size() > residues_.size());
    if (!mols_share_residues) {
      // More in-depth check
      for (std::vector<Molecule>::const_iterator mol = molecules_.begin() + 1;
                                                 mol != molecules_.end(); ++mol)
      {
        int m0_resnum = atoms_[(mol-1)->BeginAtom()].ResNum();
        int m1_resnum = atoms_[    mol->BeginAtom()].ResNum();
        if (m0_resnum == m1_resnum) {
          mols_share_residues = true;
          unsigned int molnum = mol - molecules_.begin();
          mprintf("Warning: 2 or more molecules (%u and %u) share residue numbers (%i).\n",
                  molnum, molnum+1, m0_resnum+1);
          break;
        }
      }
    }
    if (mols_share_residues) {
      //if (bondsearch)
        mprintf("Warning:   Either residue information is incorrect or molecule determination"
                " was inaccurate.\n");
      //else
      //  mprintf("Warning: Residue information appears to be incorrect.\n");
      mprintf("Warning:   Basing residue information on molecules.\n");
      std::vector<Residue> newResArray;
      unsigned int res_first_atom = 0;
      while (res_first_atom < atoms_.size()) {
        // Search for next atom with different res or molecule number.
        int current_rnum = atoms_[res_first_atom].ResNum();
        int current_mnum = atoms_[res_first_atom].MolNum();
        unsigned int res_last_atom = res_first_atom;
        while (res_last_atom != atoms_.size() &&
               atoms_[res_last_atom].ResNum() == current_rnum &&
               atoms_[res_last_atom].MolNum() == current_mnum)
          ++res_last_atom; 
        for (unsigned int r_atm = res_first_atom; r_atm != res_last_atom; ++r_atm)
          atoms_[r_atm].SetResNum( newResArray.size() ); // TODO combine with above while
        newResArray.push_back( Residue(residues_[current_rnum], res_first_atom, res_last_atom) );
        res_first_atom = res_last_atom;
      }
      mprintf("Warning:   Old # residues= %zu, new # residues = %zu\n",
              residues_.size(), newResArray.size());
      residues_ = newResArray;
      if (debug_ > 0)
        for (std::vector<Residue>::const_iterator res = newResArray.begin();
                                                  res != newResArray.end(); ++res)
          mprintf("%s first=%i last=%i orig=%i icode=%c\n",
                  res->c_str(), res->FirstAtom()+1, res->LastAtom(),
                  res->OriginalResNum(), res->Icode());

    }
  }
  // Set up solvent information
  if (SetSolventInfo())
    mprinterr("Error: Could not determine solvent information for %s.\n", c_str());
  // Determine excluded atoms
  DetermineExcludedAtoms();
  // Determine # of extra points.
  DetermineNumExtraPoints();

  return 0;
}

/** Reset any extended PDB info. */
void Topology::ResetPDBinfo() {
  int rnum = 1;
  for (std::vector<Residue>::iterator res = residues_.begin(); 
                                      res != residues_.end(); ++res, ++rnum)
  {
    res->SetOriginalNum(rnum);
    res->SetIcode(' ');
    res->SetChainID(' ');
  }
  for (std::vector<AtomExtra>::iterator ex = extra_.begin();
                                        ex != extra_.end(); ++ex)
    ex->SetAltLoc(' '); // TODO bfactor, occupancy?
}

/** For topology formats that do not contain residue info, base residues
  * on molecules.
  */
// FIXME Can the routine in CommonSetup be used in place of this instead?
int Topology::Setup_NoResInfo() {
  mprintf("\tAttempting to determine residue info from molecules.\n");
  if (DetermineMolecules()) {
    mprintf("Warning: Could not determine molecule info. Not setting up residues.\n");
    return 0;
  }
  // Save residue name if its there at all.
  NameType default_res_name, res_name;
  if (!residues_.empty())
    default_res_name = residues_[0].Name();
  else
    default_res_name = "RES";
  // Set residue info to match molecule info.
  residues_.clear();
  int resnum = 0;
  for (std::vector<Molecule>::const_iterator mol = molecules_.begin();
                                             mol != molecules_.end();
                                           ++mol, ++resnum)
  {
    // Try to detect at least water as solvent. Assume CommonSetup will be
    // run after this to set up molecule solvent info.
    if (mol->NumAtoms() == 3) {
      int nH = 0;
      int nO = 0;
      for (int atnum = mol->BeginAtom(); atnum != mol->EndAtom(); atnum++)
      {
        if (atoms_[atnum].Element() == Atom::HYDROGEN) nH++;
        if (atoms_[atnum].Element() == Atom::OXYGEN)   nO++;
      }
      if (nO == 1 && nH == 2) res_name = "HOH";
    } else
      res_name = default_res_name;
    residues_.push_back( Residue(res_name, resnum+1, ' ', ' ') );
    residues_.back().SetFirstAtom( mol->BeginAtom() );
    residues_.back().SetLastAtom( mol->EndAtom() );
    // Update atom residue numbers
    for (int atnum = residues_.back().FirstAtom(); 
             atnum != residues_.back().LastAtom(); ++atnum)
      atoms_[atnum].SetResNum( resnum );
  }
  return 0;
}

static inline int NoAtomsErr(const char* msg) {
  mprinterr("Error: Cannot set up %s, no atoms present.\n");
  return 1;
}

// Topology::Resize()
void Topology::Resize(Pointers const& pIn) {
  atoms_.clear();
  residues_.clear();
  molecules_.clear();
  radius_set_.clear();
  bonds_.clear();
  bondsh_.clear();
  bondparm_.clear();
  angles_.clear();
  anglesh_.clear();
  angleparm_.clear();
  dihedrals_.clear();
  dihedralsh_.clear();
  dihedralparm_.clear();
  nonbond_.Clear();
  cap_.Clear();
  lesparm_.Clear();
  chamber_.Clear();
  extra_.clear();
  parmBox_.SetNoBox();
  refCoords_ = Frame();
  ipol_ = 0;
  NsolventMolecules_ = 0;
  n_extra_pts_ = 0;
  n_atom_types_ = 0;

  atoms_.resize( pIn.natom_ );
  residues_.resize( pIn.nres_ );
  extra_.resize( pIn.nextra_ );
  bondparm_.resize( pIn.nBndParm_ );
  angleparm_.resize( pIn.nAngParm_ );
  dihedralparm_.resize( pIn.nDihParm_ );
}

double Topology::GetVDWradius(int a1) const {
  //TODO: return zero when no params?
  NonbondType const& LJ = GetLJparam(a1, a1);
  if (LJ.B() > 0.0)
    return ( 0.5 * pow(2.0 * LJ.A() / LJ.B(), (1.0/6.0)) );
  else
    return 0.0;
}

double Topology::GetVDWdepth(int a1) const {
  NonbondType const& LJ = GetLJparam(a1, a1);
  if (LJ.A() > 0.0)
    return ( (LJ.B() * LJ.B()) / (4.0 * LJ.A()) );
  else
    return 0.0;
}

// Topology::SetAtomBondInfo()
/** Set up bond information in the atoms array based on given BondArray.
  */
void Topology::SetAtomBondInfo(BondArray const& bonds) {
  // Add bonds based on array 
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    atoms_[ bnd->A1() ].AddBondToIdx( bnd->A2() );
    atoms_[ bnd->A2() ].AddBondToIdx( bnd->A1() );
  }
}

// -----------------------------------------------------------------------------
// Topology::AddBondParam()
/** Create parameters for given bond based on element types. */
void Topology::AddBondParam(BondType& bnd, BP_mapType& bpMap)
{
  unsigned int bp_idx;
  Atom::AtomicElementType a1Elt = atoms_[bnd.A1()].Element();
  Atom::AtomicElementType a2Elt = atoms_[bnd.A2()].Element();
  std::set<Atom::AtomicElementType> Eset;
  Eset.insert( a1Elt );
  Eset.insert( a2Elt );
  // Has this bond parameter been defined?
  BP_mapType::iterator bp = std::find(bpMap.begin(), bpMap.end(), Eset);
  if (bp == bpMap.end()) { // Bond parameter Not defined
    bp_idx = bondparm_.size();
    bpMap.push_back( Eset );
    bondparm_.push_back( BondParmType(0.0, Atom::GetBondLength(a1Elt, a2Elt)) );
  } else
    bp_idx = bp - bpMap.begin();
  //mprintf("DEBUG:\t\t%i:[%s] -- %i:[%s] Cut=%f BPidx=%u\n",
  //        bnd.A1()+1, atoms_[bnd.A1()].c_str(), bnd.A2()+1, atoms_[bnd.A2()].c_str(),
  //        bondparm_[bp_idx].Req(), bp_idx);
  bnd.SetIdx( bp_idx );
}

// Topology::AssignBondParameters()
void Topology::AssignBondParameters() {
  mprintf("Warning: %s: Determining default bond distances from element types.\n", c_str());
  bondparm_.clear();
  // Hold indices into bondparm for unique element pairs
  BP_mapType bpMap;
  for (BondArray::iterator bnd = bondsh_.begin(); bnd != bondsh_.end(); ++bnd)
    AddBondParam( *bnd, bpMap ); 
  for (BondArray::iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
    AddBondParam( *bnd, bpMap );
} 

// Topology::AddBond()
void Topology::AddBond(int atom1, int atom2, BondParmType const& BPin) {
  // See if the BondParm exists.
  int pidx = -1;
  for (BondParmArray::const_iterator bp = bondparm_.begin(); bp != bondparm_.end(); ++bp)
    if ( fabs(BPin.Rk()  - bp->Rk() ) < Constants::SMALL &&
         fabs(BPin.Req() - bp->Req()) < Constants::SMALL )
    {
      pidx = (int)(bp - bondparm_.begin());
      break;
    }
  if (pidx == -1) {
    pidx = (int)bondparm_.size();
    bondparm_.push_back( BPin );
  }
  AddBond( atom1, atom2, pidx );
}

static inline int WarnOutOfRange(int Natom, int atom, const char* type) {
  if (atom < 0 || atom >= Natom) {
    mprintf("Warning: Atom # %i is out of range, cannot create %s.\n", atom+1, type);
    return 1;
  }
  return 0;
}

// Topology::AddBond()
/** Create a bond between atom1 and atom2, update the atoms array.
  * For bonds to H always insert the H second.
  */
void Topology::AddBond(int atom1, int atom2, int pidxIn) {
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), atom1, "bond")) return;
  if (WarnOutOfRange(atoms_.size(), atom2, "bond")) return;
  // Check for duplicate bond
  for (Atom::bond_iterator ba = atoms_[atom1].bondbegin();
                           ba != atoms_[atom1].bondend(); ++ba)
    if ( *ba == atom2 ) {
      if (debug_ > 0)
        mprintf("Warning: Bond between atoms %i and %i already exists.\n", atom1+1, atom2+1);
      return;
    }
  // Check if parm index is out of range;
  int pidx;
  if (pidxIn < (int)bondparm_.size())
    pidx = pidxIn;
  else {
    mprintf("Warning: No bond parameters for index %i\n", pidxIn);
    pidx = -1;
  }
  bool a1H = (atoms_[atom1].Element() == Atom::HYDROGEN);
  bool a2H = (atoms_[atom2].Element() == Atom::HYDROGEN);
  //mprintf("\t\t\tAdding bond %i to %i (isH=%i)\n",atom1+1,atom2+1,(int)isH);
  // Update bonds arrays
  if (a1H || a2H) {
    if (a1H)
      bondsh_.push_back( BondType(atom2, atom1, pidx) );
    else
      bondsh_.push_back( BondType(atom1, atom2, pidx) );
  } else
    bonds_.push_back( BondType( atom1, atom2, pidx ) );
  // Update atoms
  atoms_[atom1].AddBondToIdx( atom2 );
  atoms_[atom2].AddBondToIdx( atom1 );
}

/** For use when element data may not yet be available. If isH, it is
  * assumed that the second atom is the H.
  */
void Topology::AddBond(BondType const& bndIn, bool isH) {
  if (isH)
    bondsh_.push_back( bndIn );
  else
    bonds_.push_back( bndIn );
  // Update atoms
  atoms_[bndIn.A1()].AddBondToIdx( bndIn.A2() );
  atoms_[bndIn.A2()].AddBondToIdx( bndIn.A1() );
}

// Topology::AddAngle() 
void Topology::AddAngle(int atom1, int atom2, int atom3, AngleParmType const& APin) {
  // See if the AngleParm exists.
  int pidx = -1;
  for (AngleParmArray::const_iterator ap = angleparm_.begin(); ap != angleparm_.end(); ++ap)
    if ( fabs(APin.Tk()  - ap->Tk() ) < Constants::SMALL &&
         fabs(APin.Teq() - ap->Teq()) < Constants::SMALL )
    {
      pidx = (int)(ap - angleparm_.begin());
      break;
    }
  if (pidx == -1) {
    pidx = (int)angleparm_.size();
    angleparm_.push_back( APin );
  }
  AddAngle( atom1, atom2, atom3, pidx );
}

// Topology::AddAngle() 
void Topology::AddAngle(int atom1, int atom2, int atom3, int pidxIn) {
  // FIXME: Check duplicate
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), atom1, "angle")) return;
  if (WarnOutOfRange(atoms_.size(), atom2, "angle")) return;
  if (WarnOutOfRange(atoms_.size(), atom3, "angle")) return;
  // Check if parm index is out of range;
  int pidx;
  if (pidxIn < (int)angleparm_.size())
    pidx = pidxIn;
  else {
    mprintf("Warning: No angle parameters for index %i\n", pidxIn);
    pidx = -1;
  }
  // Update angle arrays
  if (atoms_[atom1].Element() == Atom::HYDROGEN ||
      atoms_[atom2].Element() == Atom::HYDROGEN ||
      atoms_[atom3].Element() == Atom::HYDROGEN)
    anglesh_.push_back( AngleType(atom1, atom2, atom3, pidx) );
  else
    angles_.push_back( AngleType(atom1, atom2, atom3, pidx) );
}

// Topology::AddAngle()
void Topology::AddAngle(AngleType const& angIn, bool isH) {
  if (isH)
    anglesh_.push_back( angIn );
  else
    angles_.push_back( angIn );
}

// Topology::AddDihedral()
void Topology::AddDihedral(DihedralType const& dih, DihedralParmType const& DPin)
{
  // See if the DihedralParm exists.
  int pidx = -1;
  for (DihedralParmArray::const_iterator dp = dihedralparm_.begin();
                                         dp != dihedralparm_.end(); ++dp)
    if ( fabs(DPin.Pk()    - dp->Pk()   ) < Constants::SMALL &&
         fabs(DPin.Pn()    - dp->Pn()   ) < Constants::SMALL &&
         fabs(DPin.Phase() - dp->Phase()) < Constants::SMALL &&
         fabs(DPin.SCEE()  - dp->SCEE() ) < Constants::SMALL &&
         fabs(DPin.SCNB()  - dp->SCNB() ) < Constants::SMALL )
    {
      pidx = (int)(dp - dihedralparm_.begin());
      break;
    }
  if (pidx == -1) {
    pidx = (int)dihedralparm_.size();
    dihedralparm_.push_back( DPin );
  }
  AddDihedral( dih, pidx );
}

// Topology::AddDihedral()
void Topology::AddDihedral(DihedralType const& dihIn, int pidxIn) {
  // FIXME: Check duplicate
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), dihIn.A1(), "dihedral")) return;
  if (WarnOutOfRange(atoms_.size(), dihIn.A2(), "dihedral")) return;
  if (WarnOutOfRange(atoms_.size(), dihIn.A3(), "dihedral")) return;
  if (WarnOutOfRange(atoms_.size(), dihIn.A4(), "dihedral")) return;
  // Check if parm index is out of range;
  int pidx;
  if (pidxIn < (int)dihedralparm_.size())
    pidx = pidxIn;
  else {
    mprintf("Warning: No dihedral parameters for index %i\n", pidxIn);
    pidx = -1;
  }
  DihedralType dih = dihIn;
  dih.SetIdx( pidx );
  // Update dihedral arrays
  if (atoms_[dih.A1()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A2()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A3()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A4()].Element() == Atom::HYDROGEN)
    dihedralsh_.push_back( dih );
  else
    dihedrals_.push_back( dih );
}

void Topology::AddDihedral(DihedralType const& dihIn, bool isH) {
  if (isH)
    dihedralsh_.push_back( dihIn );
  else
    dihedrals_.push_back( dihIn );
}

// -----------------------------------------------------------------------------
// Topology::VisitAtom()
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

/** Recursive search for molecules along bonds of each atom. */
int Topology::RecursiveMolSearch() {
  //Timer t_stack;
  //t_stack.Start();
  int atomnum = 0;
  int mol = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin(); atom != atoms_.end(); atom++)
  {
    if ( atom->NoMol() ) {
      VisitAtom( atomnum, mol );
      ++mol;
    }
    ++atomnum;
  }
  //t_stack.Stop();
  //t_stack.WriteTiming(1, "Recursive mol search:");
  return mol;
}

/** Non-recursive molecule search. Better for larger systems, uses the heap. */
int Topology::NonrecursiveMolSearch() {
  if (debug_ > 0) mprintf("DEBUG: Beginning non-recursive molecule search.\n");
  // Recursive search for high atom counts can blow the stack away.
  //Timer t_nostack;
  //t_nostack.Start();
  std::stack<unsigned int> nextAtomToSearch;
  bool unassignedAtomsRemain = true;
  unsigned int currentAtom = 0;
  unsigned int currentMol = 0;
  unsigned int lowestUnassignedAtom = 0;
  while (unassignedAtomsRemain) {
    // This atom is in molecule.
    atoms_[currentAtom].SetMol( currentMol );
    //mprintf("DEBUG:\tAssigned atom %u to mol %u\n", currentAtom, currentMol);
    // All atoms bonded to this one are in molecule.
    for (Atom::bond_iterator batom = atoms_[currentAtom].bondbegin();
                             batom != atoms_[currentAtom].bondend(); ++batom)
    {
      if (atoms_[*batom].NoMol()) {
        if (atoms_[*batom].Nbonds() > 1)
          // Bonded atom has more than 1 bond; needs to be searched.
          nextAtomToSearch.push( *batom );
        else {
          // Bonded atom only bonded to current atom. No more search needed.
          atoms_[*batom].SetMol( currentMol );
          //mprintf("DEBUG:\t\tAssigned bonded atom %i to mol %u\n", *batom, currentMol);
        }
      }
    }
    if (nextAtomToSearch.empty()) {
      //mprintf("DEBUG:\tNo atoms left in stack. Searching for next unmarked atom.\n");
      // No more atoms to search. Find next unmarked atom.
      currentMol++;
      unsigned int idx = lowestUnassignedAtom;
      for (; idx != atoms_.size(); idx++)
        if (atoms_[idx].NoMol()) break;
      if (idx == atoms_.size())
        unassignedAtomsRemain = false;
      else {
        currentAtom = idx;
        lowestUnassignedAtom = idx + 1;
      }
    } else {
      currentAtom = nextAtomToSearch.top();
      nextAtomToSearch.pop();
      //mprintf("DEBUG:\tNext atom from stack: %u\n", currentAtom);
    }
  }
  //t_nostack.Stop();
  //t_nostack.WriteTiming(1, "Non-recursive mol search:");
  return (int)currentMol;
}

// Topology::ClearMolecules()
/** Clear molecules and reset molecule info for each atom. */
void Topology::ClearMolecules() {
  molecules_.clear();
  for (std::vector<Atom>::iterator atom = atoms_.begin(); atom != atoms_.end(); atom++)
    atom->SetMol( -1 );
}

// Topology::DetermineMolecules()
/** Determine individual molecules using bond information. Performs a 
  * recursive search over the bonds of each atom.
  */
int Topology::DetermineMolecules() {
  // Since this is always done only print when debugging
  if (debug_>0) mprintf("\t%s: determining molecule info from bonds.\n",c_str());
  // Reset molecule info for each atom
  ClearMolecules();
  int numberOfMolecules = 0;
  if (atoms_.size() > 150000) // Seems to be when performance of nonrecursive approaches recursive
    numberOfMolecules = NonrecursiveMolSearch();
  else
    numberOfMolecules = RecursiveMolSearch();
/*// DEBUG Compare both methods
  int test_nmol = NonrecursiveMolSearch();
  std::vector<int> molNums( atoms_.size() );
  for (unsigned int idx = 0; idx != atoms_.size(); idx++)
    molNums[idx] = atoms_[idx].MolNum();
  ClearMolecules();
  numberOfMolecules = RecursiveMolSearch();
  if (test_nmol != numberOfMolecules)
    mprintf("Num mols found with non-recursive search (%i) does not match (%i)\n",
            test_nmol, numberOfMolecules);
  for (unsigned int idx = 0; idx != atoms_.size(); idx++)
    if (molNums[idx] != atoms_[idx].MolNum())
      mprintf("%u: Mol num in non-recursive search %i does not match %i\n",
              idx, molNums[idx], atoms_[idx].MolNum());
*/
  if (debug_ > 0) {
    mprintf("\t%i molecules.\n", numberOfMolecules);
    if (debug_ > 1)
    for (std::vector<Atom>::const_iterator atom = atoms_.begin(); atom != atoms_.end(); ++atom)
      mprintf("\t\tAtom %i assigned to molecule %i\n", atom - atoms_.begin(), atom->MolNum());
  }

  // Update molecule information
  molecules_.resize( numberOfMolecules );
  if (numberOfMolecules == 0) return 0;
  std::vector<Molecule>::iterator molecule = molecules_.begin();
  molecule->SetFirst(0);
  std::vector<Atom>::const_iterator atom = atoms_.begin(); 
  int lastMol = atom->MolNum();
  int atomNum = 0;
  for (; atom != atoms_.end(); atom++)
  {
    if ( atom->MolNum() > lastMol ) {
      // Set last atom of molecule
      molecule->SetLast( atomNum );
      // Set first atom of next molecule
      ++molecule;
      molecule->SetFirst( atomNum );
      lastMol = atom->MolNum();
    } else if ( atom->MolNum()  < lastMol) {
      mprinterr("Error: Atom %u was assigned a lower molecule # than previous atom.\n"
                "Error:   This can happen if bond information is incorrect or missing, or if the\n"
                "Error:   atom numbering in molecules is not sequential. Try one of the\n"
                "Error:   following:\n"
                "Error: - If this is a PDB file, try using the 'noconect' keyword.\n"
                "Error: - If this topology did not have bond info, try increasing the bond\n"
                "Error:   search cutoff above 0.2 Ang. ('bondsearch <cutoff>').\n"
                "Error: - Use the 'fixatomorder' command to reorder the topology and any\n"
                "Error:   associated coordinates.\n"
                "Error: - Use the 'setMolecules' command in parmed to reorder only the\n"
                "Error:   topology.\n", atom - atoms_.begin() + 1);
      ClearMolecules();
      return 1;
    }
    ++atomNum;
  }
  molecule->SetLast( atoms_.size() );
  return 0;
}

// -----------------------------------------------------------------------------
// Topology::AtomDistance()
void Topology::AtomDistance(int originalAtom, int atom, int dist, std::set<int> &excluded) const 
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
    //mprintf("    Determining excluded atoms for atom %i\n",atomi+1);
    // AtomDistance recursively sets each atom bond distance from atomi
    AtomDistance(atomi, atomi, 0, excluded_i);
    atoms_[atomi].AddExclusionList( excluded_i );
    // DEBUG
    //mprintf("\tAtom %i Excluded:",atomi+1);
    //for (Atom::excluded_iterator ei = atoms_[atomi].excludedbegin(); 
    //                             ei != atoms_[atomi].excludedend(); ++ei)
    //  mprintf(" %i",*ei + 1);
    //mprintf("\n");
  } // END loop over atomi
}

// Topology::DetermineNumExtraPoints()
void Topology::DetermineNumExtraPoints() {
  n_extra_pts_ = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
    if ( (*atom).Element() == Atom::EXTRAPT ) ++n_extra_pts_;
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
  // If maskexpr is empty this means remove all solvent information.
  if (maskexpr.empty()) {
    mprintf("Warning: Removing all solvent information from %s\n", c_str());
    for (std::vector<Molecule>::iterator mol = molecules_.begin(); 
                                         mol != molecules_.end(); ++mol)
      mol->SetNoSolvent();
    NsolventMolecules_ = 0;
    return 0;
  }
  // Setup mask
  CharMask mask( maskexpr );
  SetupCharMask( mask );
  if (mask.None()) {
    mprinterr("Error: SetSolvent [%s]: Mask %s selects no atoms.\n", c_str(), maskexpr.c_str());
    return 1;
  }
  // Loop over all molecules
  NsolventMolecules_ = 0;
  int numSolvAtoms = 0;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); ++mol)
  {
    // Reset old solvent information.
    mol->SetNoSolvent();
    // If any atoms in this molecule are selected by mask, make entire
    // molecule solvent.
    for (int atom = mol->BeginAtom(); atom < mol->EndAtom(); ++atom) {
      if ( mask.AtomInCharMask( atom ) ) {
        mol->SetSolvent();
        ++NsolventMolecules_;
        numSolvAtoms += mol->NumAtoms();
        break;
      }
    }
  }

  mprintf("\tSolvent Mask [%s]: %i solvent molecules, %i solvent atoms\n",
          maskexpr.c_str(), NsolventMolecules_, numSolvAtoms);
  return 0;
}

// Topology::SetSolventInfo()
/** Determine which molecules are solvent based on residue name. */
int Topology::SetSolventInfo() {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolventInfo: No molecule information.\n");
    return 1;
  }
  // Loop over each molecule. Check if first residue of molecule is solvent.
  NsolventMolecules_ = 0;
  int numSolvAtoms = 0;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); mol++)
  {
    int firstRes = atoms_[ mol->BeginAtom() ].ResNum();
    if ( residues_[firstRes].NameIsSolvent() ) {
      mol->SetSolvent();
      ++NsolventMolecules_;
      numSolvAtoms += mol->NumAtoms();
    }
  }

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
int Topology::SetupIntegerMask(AtomMask &mask) const {
  return mask.SetupMask(atoms_, residues_, molecules_, refCoords_.xAddress());
}

// Topology::SetupCharMask()
int Topology::SetupCharMask(CharMask &mask) const {
  return mask.SetupMask(atoms_, residues_, molecules_, refCoords_.xAddress());
}

// Topology::SetupIntegerMask()
int Topology::SetupIntegerMask(AtomMask &mask, Frame const& frame) const {
  if (frame.empty()) return mask.SetupMask(atoms_, residues_, molecules_, 0);
  return mask.SetupMask(atoms_, residues_, molecules_, frame.xAddress());
}

// Topology::SetupCharMask()
int Topology::SetupCharMask(CharMask &mask, Frame const& frame) const {
  if (frame.empty()) return mask.SetupMask(atoms_, residues_, molecules_, 0);
  return mask.SetupMask(atoms_, residues_, molecules_, frame.xAddress());
}

// -----------------------------------------------------------------------------
int Topology::scale_dihedral_K(DihedralArray& dihedrals, CharMask const& Mask,
                               double scale_factor, bool useAll)
{
  std::vector<int> newDihedralParms( dihedralparm_.size(), -1 );
  for (DihedralArray::iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    bool validDihedral;
    if (useAll)
      validDihedral= ( Mask.AtomInCharMask(dih->A1()) && Mask.AtomInCharMask(dih->A2()) &&
                       Mask.AtomInCharMask(dih->A3()) && Mask.AtomInCharMask(dih->A4()) );
    else
      validDihedral= ( Mask.AtomInCharMask(dih->A1()) || Mask.AtomInCharMask(dih->A2()) || 
                       Mask.AtomInCharMask(dih->A3()) || Mask.AtomInCharMask(dih->A4()) );
    if (validDihedral) {
      // See if this dihedral type was previously scaled.
      int oldidx = dih->Idx();
      if (oldidx == -1) {
        mprinterr("Error: No dihedral parameters.\n");
        return 1;
      }
      int newidx = newDihedralParms[oldidx];
      if (newidx == -1) {
        // Scale and add new dihedral parameter type.
        DihedralParmType newparm = dihedralparm_[oldidx];
        newparm.Pk() *= scale_factor;
        newidx = (int)dihedralparm_.size();
        dihedralparm_.push_back( newparm );
        newDihedralParms[oldidx] = newidx;
      } 
      // Update dihedral parameter index.
      dih->SetIdx( newidx );
      mprintf("\tDihedral %s-%s-%s-%s old PK= %g  new PK= %g\n",
              AtomMaskName(dih->A1()).c_str(),
              AtomMaskName(dih->A2()).c_str(),
              AtomMaskName(dih->A3()).c_str(),
              AtomMaskName(dih->A4()).c_str(),
              dihedralparm_[oldidx].Pk(), dihedralparm_[newidx].Pk());
    }
  }
  return 0;
}

int Topology::ScaleDihedralK(double scale_factor, std::string const& maskExpr, bool useAll)
{
  if (maskExpr.empty()) {
    // Scale all
    for (DihedralParmArray::iterator dk = dihedralparm_.begin();
                                     dk != dihedralparm_.end(); ++dk)
      dk->Pk() *= scale_factor;
  } else {
    // Scale only dihedrals with atoms in mask. Requires adding new types.
    CharMask Mask( maskExpr );
    if (SetupCharMask( Mask )) return 1;
    if (scale_dihedral_K( dihedrals_,  Mask, scale_factor, useAll )) return 1;
    if (scale_dihedral_K( dihedralsh_, Mask, scale_factor, useAll )) return 1;
  }
  return 0;
}

// Topology::ModifyByMap()
/** \return Pointer to new Topology based on this Topology, deleting atoms
  *         that are not in the given map (Map[newatom] = oldatom).
  */
Topology* Topology::ModifyByMap(std::vector<int> const& MapIn, bool setupFullParm) const {
  Topology *newParm = new Topology();

  newParm->parmName_ = parmName_;
  newParm->fileName_ = fileName_;
  newParm->radius_set_ = radius_set_;
  newParm->debug_ = debug_;
  newParm->n_atom_types_ = n_atom_types_;

  // Reverse Atom map
  std::vector<int> atomMap( atoms_.size(),-1 );
  // Save solvent status of atoms
  std::vector<bool> isSolvent;
  isSolvent.reserve( MapIn.size() );

  // Copy atoms from this parm that are in Mask to newParm.
  int oldres = -1;
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
      Residue const& cr = residues_[curres];
      newParm->residues_.push_back( Residue(cr.Name(), cr.OriginalResNum(),
                                            cr.Icode(), cr.ChainID()) );
      newParm->residues_.back().SetFirstAtom( newatom );
      oldres = curres;
    }
    // Clear bond information from new atom
    newparmAtom.ClearBonds();
    // Set new atom num and residue num
    newparmAtom.SetResNum( newParm->residues_.size() - 1 );
    // Check if this atom belongs to a solvent molecule.
    if (!molecules_.empty())
      isSolvent.push_back( Mol(newparmAtom.MolNum()).IsSolvent() );
    // Place new atom in newParm
    newParm->atoms_.push_back( newparmAtom );
  }
  if (newParm->atoms_.empty()) {
    mprintf("Warning: All atoms have been stripped.\n");
    return newParm;
  }

  // Set last residue last atom
  newParm->residues_.back().SetLastAtom( newParm->atoms_.size() );

  // Copy reference if present
  if (!refCoords_.empty()) {
    newParm->refCoords_.SetupFrameM( atoms_ );
    newParm->refCoords_.ModifyByMap( refCoords_, MapIn );
  }

  // NOTE: Since in the bond/angle/dihedral atom arrays the parm indices have 
  //       survived intact we can just include direct copies of all the 
  //       parameter arrays for now. May want to cull unused params later.

  // Set up new bond information
  newParm->bonds_ = StripBondArray( bonds_, atomMap );
  newParm->bondsh_ = StripBondArray( bondsh_, atomMap );
  newParm->SetAtomBondInfo( newParm->bonds_ );
  newParm->SetAtomBondInfo( newParm->bondsh_ );
  std::vector<int> parmMap( bondparm_.size(), -1 ); // Map[oldidx] = newidx
  StripBondParmArray( newParm->bonds_,  parmMap, newParm->bondparm_ );
  StripBondParmArray( newParm->bondsh_, parmMap, newParm->bondparm_ );
  //mprintf("DEBUG: Original bond parm array= %zu, new bond parm array = %zu\n",
  //        bondparm_.size(), newParm->bondparm_.size());
  // Give stripped parm the same pindex as original
  newParm->pindex_ = pindex_;
  // Copy box information
  newParm->parmBox_ = parmBox_;
  // If we dont care about setting up full parm information, exit now.
  if (!setupFullParm) return newParm;

  // Set new molecule information based on new bonds
  if (newParm->DetermineMolecules()) {
    mprintf("Warning: Could not set up molecule information for stripped topology %s\n",
            newParm->c_str());
  }

  // Determine solvent
  if (!molecules_.empty()) {
    // Set new solvent information based on old molecules
    // For speed just check the first atom. A strip should never create new
    // molecules, just break up previous ones, so if a new molecule has an atom
    // that was solvent it is still solvent.
    newParm->NsolventMolecules_ = 0;
    for (std::vector<Molecule>::iterator mol = newParm->molecules_.begin();
                                         mol != newParm->molecules_.end(); ++mol)
    {
      if ( isSolvent[ mol->BeginAtom() ] ) {
        mol->SetSolvent();
        newParm->NsolventMolecules_++;
      } else
        mol->SetNoSolvent();
    }
  } else {
    // No solvent information previously. Check if there is solvent now,
    // which could be the case if molecule information was not previously
    // determined.
    newParm->SetSolventInfo();
  }

  // Set up new angle info
  newParm->angles_ = StripAngleArray( angles_, atomMap );
  newParm->anglesh_ = StripAngleArray( anglesh_, atomMap );
  if (!angleparm_.empty()) {
    parmMap.assign( angleparm_.size(), -1 );
    StripAngleParmArray( newParm->angles_,  parmMap, newParm->angleparm_ );
    StripAngleParmArray( newParm->anglesh_, parmMap, newParm->angleparm_ );
  }
  // Set up new dihedral info
  newParm->dihedrals_ = StripDihedralArray( dihedrals_, atomMap );
  newParm->dihedralsh_ = StripDihedralArray( dihedralsh_, atomMap );
  if (!dihedralparm_.empty()) {
    parmMap.assign( dihedralparm_.size(), -1 );
    StripDihedralParmArray( newParm->dihedrals_,  parmMap, newParm->dihedralparm_ );
    StripDihedralParmArray( newParm->dihedralsh_, parmMap, newParm->dihedralparm_ );
  }
  // Set up nonbond info. First determine which atom types remain.
  if (nonbond_.HasNonbond()) {
    parmMap.clear();               // parmMap[oldtype]      = newtype
    std::vector<int> oldTypeArray; // oldTypeArray[newtype] = oldtype
    for (std::vector<Atom>::const_iterator atm = newParm->atoms_.begin();
                                           atm != newParm->atoms_.end(); ++atm)
    {
      int oldidx = atm->TypeIndex();
      if (oldidx >= (int)parmMap.size())
        parmMap.resize( oldidx+1, -1 );
      if (parmMap[oldidx] == -1) {
        parmMap[oldidx] = (int)oldTypeArray.size();
        oldTypeArray.push_back( oldidx );
      }
      //int newidx = parmMap[oldidx];
      //mprintf("DEBUG: '%s' Old type index=%i, new type index = %i\n", atm->c_str(), oldidx, newidx);
    }
    //mprintf("DEBUG: # new types %zu\n", oldTypeArray.size());
    // Set up new nonbond and nonbond index arrays.
    newParm->nonbond_.SetNtypes( oldTypeArray.size() );
    if (chamber_.HasChamber())
      newParm->chamber_.SetNLJ14terms( (oldTypeArray.size()*(oldTypeArray.size()+1))/2 );
    for (int a1idx = 0; a1idx != (int)oldTypeArray.size(); a1idx++)
    {
      int atm1 = oldTypeArray[a1idx];
      for (int a2idx = a1idx; a2idx != (int)oldTypeArray.size(); a2idx++)
      {
        int atm2 = oldTypeArray[a2idx];
        int oldnbidx = nonbond_.GetLJindex( atm1, atm2 );
        if (oldnbidx > -1) {
          // This is a traditional LJ 6-12 term. Because of the way the LJ 1-4
          // code is laid out in sander/pmemd the LJ matrix has to be laid out
          // indepdendent of the nonbond index array.
          newParm->nonbond_.AddLJterm( a1idx, a2idx, nonbond_.NBarray(oldnbidx) );
        } else {
          // This is an old LJ 10-12 hbond term. Add one to the LJ 6-12 matrix
          // and one to the hbond since that seems to be the convention.
          newParm->nonbond_.AddLJterm( a1idx, a2idx, NonbondType() );
          newParm->nonbond_.AddHBterm( a1idx, a2idx, nonbond_.HBarray((-oldnbidx)-1) );
        }
        //int newnbidx = newParm->nonbond_.GetLJindex( a1idx, a2idx );
        //mprintf("DEBUG: oldtypei=%i oldtypej=%i Old NB index=%i, newtypi=%i newtypej=%i new NB idx=%i testidx=%i\n", 
        //        atm1, atm2, oldnbidx, a1idx, a2idx, newnbidx, testidx);
        if (chamber_.HasChamber()) {
          // Update LJ 1-4 as well. No need to worry about hbond terms here,
          // just recalculate the old index and determine new one.
          int ibig = std::max(atm1, atm2) + 1;
          int isml = std::min(atm1, atm2) + 1;
              oldnbidx = (ibig*(ibig-1)/2+isml)-1;
              ibig = a2idx + 1;
              isml = a1idx + 1;
          int newnbidx = (ibig*(ibig-1)/2+isml)-1;
          newParm->chamber_.SetLJ14( newnbidx ) = chamber_.LJ14()[oldnbidx];
        }
      }
    }
    // Update atom type indices.
    for (std::vector<Atom>::iterator atm = newParm->atoms_.begin();
                                     atm != newParm->atoms_.end(); ++atm)
      atm->SetTypeIndex( parmMap[atm->TypeIndex()] );
  }
  // LES info - FIXME: Not sure if stripping this is valid so print a warning.
  if (lesparm_.HasLES()) {
    mprintf("Warning: LES info present. Stripped topology may not have correct LES info.\n");
    newParm->lesparm_.SetTypes( lesparm_.Ntypes(), lesparm_.FAC() );
    for (std::vector<int>::const_iterator old_it = MapIn.begin(); old_it != MapIn.end(); ++old_it)
    {
      if (*old_it >= 0)
        newParm->lesparm_.AddLES_Atom( lesparm_.Array()[*old_it] );
    }
  }
  // CAP info - dont support stripping such topologies right now
  if (cap_.HasWaterCap())
    mprintf("Warning: Stripping of CAP info not supported. Removing CAP info.\n");
  // CHAMBER info
  if (chamber_.HasChamber()) {
    newParm->chamber_.SetHasChamber( true );
    newParm->chamber_.SetDescription( chamber_.Description() );
    // Urey-Bradley
    newParm->chamber_.SetUB() = StripBondArray(chamber_.UB(),atomMap);
    parmMap.assign( chamber_.UBparm().size(), -1 ); // Map[oldidx] = newidx
    StripBondParmArray( newParm->chamber_.SetUB(), parmMap,
                        newParm->chamber_.SetUBparm(), chamber_.UBparm() );
    // Impropers
    newParm->chamber_.SetImpropers() = StripDihedralArray(chamber_.Impropers(), atomMap);
    parmMap.assign( chamber_.ImproperParm().size(), -1 );
    StripDihedralParmArray( newParm->chamber_.SetImpropers(), parmMap,
                            newParm->chamber_.SetImproperParm(), chamber_.ImproperParm() );
    // NOTE 1-4 LJ parameters handled above
    // CMAP terms
    if (chamber_.HasCmap()) {
      // NOTE that atom indexing is updated but cmap indexing is not. So if
      // any CMAP terms remain all CMAP entries remain.
      for (CmapArray::const_iterator cmap = chamber_.Cmap().begin();
                                     cmap != chamber_.Cmap().end(); ++cmap)
      {
        int newA1 = atomMap[ cmap->A1() ];
        if (newA1 != -1) {
          int newA2 = atomMap[ cmap->A2() ];
          if (newA2 != -1) {
            int newA3 = atomMap[ cmap->A3() ];
            if (newA3 != -1) {
              int newA4 = atomMap[ cmap->A4() ];
              if (newA4 != -1) {
                int newA5 = atomMap[ cmap->A5() ];
                if (newA5 != -1)
                  newParm->chamber_.AddCmapTerm( CmapType(newA1,newA2,newA3,
                                                          newA4,newA5,cmap->Idx()) );
              }
            }
          }
        }
      }
      // Only add CMAP grids if there are CMAP terms left.
      if (!newParm->chamber_.Cmap().empty()) {
        for (CmapGridArray::const_iterator g = chamber_.CmapGrid().begin();
                                           g != chamber_.CmapGrid().end(); ++g)
          newParm->chamber_.AddCmapGrid( *g );
      }
    }
  }
  // Amber extra info.
  if (!extra_.empty()) {
    for (std::vector<int>::const_iterator old_it = MapIn.begin(); old_it != MapIn.end(); ++old_it)
      if (*old_it >= 0)
        newParm->extra_.push_back( extra_[*old_it] );
  }
  
  // Setup excluded atoms list - Necessary?
  newParm->DetermineExcludedAtoms();

  // Determine number of extra points
  newParm->DetermineNumExtraPoints();

  return newParm;
}

/** \return BondArray with bonds for which both atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
BondArray Topology::StripBondArray(BondArray const& bondsIn, std::vector<int> const& atomMap) const {
  BondArray bondsOut;
  // Go through old array. Use atomMap to determine what goes into newArray.
  for (BondArray::const_iterator oldbond = bondsIn.begin(); oldbond != bondsIn.end(); ++oldbond) {
    int newA1 = atomMap[ oldbond->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ oldbond->A2() ];
      if (newA2 != -1)
        bondsOut.push_back( BondType(newA1, newA2, oldbond->Idx() ) );
    }
  }
  return bondsOut;
}

/** \return AngleArray with angles for which all atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
AngleArray Topology::StripAngleArray(AngleArray const& anglesIn, std::vector<int> const& atomMap) const {
  AngleArray anglesOut;
  for (AngleArray::const_iterator oldangle = anglesIn.begin(); oldangle != anglesIn.end(); ++oldangle) {
    int newA1 = atomMap[ oldangle->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ oldangle->A2() ];
      if (newA2 != -1) {
        int newA3 = atomMap[ oldangle->A3() ];
        if (newA3 != -1)
          anglesOut.push_back( AngleType(newA1, newA2, newA3, oldangle->Idx()) );
      }
    }
  }
  return anglesOut;
}

/** \return DihedralArray with dihedrals for which all atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
DihedralArray Topology::StripDihedralArray(DihedralArray const& dihIn, std::vector<int> const& atomMap) const {
  DihedralArray dihOut;
  for (DihedralArray::const_iterator olddih = dihIn.begin(); olddih != dihIn.end(); ++olddih) {
    int newA1 = atomMap[ olddih->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ olddih->A2() ];
      if (newA2 != -1) {
        int newA3 = atomMap[ olddih->A3() ];
        if (newA3 != -1) {
          int newA4 = atomMap[ olddih->A4() ];
          if (newA4 != -1) {
            // Since in Amber improper/end dihedrals are stored as negative #s,
            // atom index 0 cannot be in 3rd or 4th position. Reverse.
            if (olddih->Type() != DihedralType::NORMAL && (newA3 == 0 || newA4 == 0))
              dihOut.push_back( DihedralType( newA4, newA3, newA2, newA1, 
                                              olddih->Type(), olddih->Idx() ) );
            else
              dihOut.push_back( DihedralType( newA1, newA2, newA3, newA4, 
                                              olddih->Type(), olddih->Idx() ) );
          }
        }
      }
    }
  }
  return dihOut;
}

// Topology::StripBondParmArray()
void Topology::StripBondParmArray(BondArray& newBondArray, std::vector<int>& parmMap,
                                  BondParmArray& newBondParm) const
{
  StripBondParmArray(newBondArray, parmMap, newBondParm, bondparm_);
}

// Topology::StripBondParmArray()
void Topology::StripBondParmArray(BondArray& newBondArray, std::vector<int>& parmMap,
                                  BondParmArray& newBondParm,
                                  BondParmArray const& oldParm) const
{
  for (BondArray::iterator bnd = newBondArray.begin();
                           bnd != newBondArray.end(); ++bnd)
  {
    int oldidx = bnd->Idx();
    int newidx = parmMap[bnd->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newBondParm.size();
      parmMap[oldidx] = newidx;
      newBondParm.push_back( oldParm[oldidx] );
    }
    //mprintf("DEBUG: Old bond parm index=%i, new bond parm index=%i\n", oldidx, newidx);
    bnd->SetIdx( newidx );
  }
}

// Topology::StripAngleParmArray()
void Topology::StripAngleParmArray(AngleArray& newAngleArray, std::vector<int>& parmMap,
                                   AngleParmArray& newAngleParm) const
{
  for (AngleArray::iterator ang = newAngleArray.begin();
                            ang != newAngleArray.end(); ++ang)
  {
    int oldidx = ang->Idx();
    int newidx = parmMap[ang->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newAngleParm.size();
      parmMap[oldidx] = newidx;
      newAngleParm.push_back( angleparm_[oldidx] );
    }
    //mprintf("DEBUG: Old angle parm index=%i, new angle parm index=%i\n", oldidx, newidx);
    ang->SetIdx( newidx );
  }
}

// Topology::StripDihedralParmArray()
void Topology::StripDihedralParmArray(DihedralArray& newDihedralArray, std::vector<int>& parmMap,
                                      DihedralParmArray& newDihedralParm) const
{
  StripDihedralParmArray(newDihedralArray, parmMap, newDihedralParm, dihedralparm_);
}

// Topology::StripDihedralParmArray()
/** Create new dihedral parm array from old one; update indices in dihedral array. */
void Topology::StripDihedralParmArray(DihedralArray& newDihedralArray, std::vector<int>& parmMap,
                                      DihedralParmArray& newDihedralParm,
                                      DihedralParmArray const& oldParm) const
{
  for (DihedralArray::iterator dih = newDihedralArray.begin();
                               dih != newDihedralArray.end(); ++dih)
  {
    int oldidx = dih->Idx();
    int newidx = parmMap[dih->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newDihedralParm.size();
      parmMap[oldidx] = newidx;
      newDihedralParm.push_back( oldParm[oldidx] );
    }
    //mprintf("DEBUG: Old dihedral parm index=%i, new dihedral parm index=%i\n", oldidx, newidx);
    dih->SetIdx( newidx );
  }
}

// Topology::AddBondArray()
void Topology::AddBondArray(BondArray const& barray, BondParmArray const& bp, int atomOffset) {
  if (bp.empty())
    for (BondArray::const_iterator bond = barray.begin(); bond != barray.end(); ++bond)
      AddBond( bond->A1() + atomOffset, bond->A2() + atomOffset );
  else
    for (BondArray::const_iterator bond = barray.begin(); bond != barray.end(); ++bond)
      AddBond( bond->A1() + atomOffset, bond->A2() + atomOffset, bp[bond->Idx()] );
}

// Topology::AddAngleArray()
void Topology::AddAngleArray(AngleArray const& aarray, AngleParmArray const& ap, int atomOffset) {
  if (ap.empty())
    for (AngleArray::const_iterator angle = aarray.begin(); angle != aarray.end(); ++angle)
      AddAngle( angle->A1() + atomOffset,
                angle->A2() + atomOffset,
                angle->A3() + atomOffset );
  else
    for (AngleArray::const_iterator angle = aarray.begin(); angle != aarray.end(); ++angle)
      AddAngle( angle->A1() + atomOffset,
                angle->A2() + atomOffset,
                angle->A3() + atomOffset, ap[angle->Idx()] );
}

// Topology::AddDihArray()
void Topology::AddDihArray(DihedralArray const& darray, DihedralParmArray const& dp, int atomOffset)
{
  if (dp.empty())
    for (DihedralArray::const_iterator dih = darray.begin(); dih != darray.end(); ++dih)
      AddDihedral( DihedralType( dih->A1() + atomOffset, dih->A2() + atomOffset,
                                 dih->A3() + atomOffset, dih->A4() + atomOffset,
                                 dih->Type() ), -1 );
  else
    for (DihedralArray::const_iterator dih = darray.begin(); dih != darray.end(); ++dih)
      AddDihedral( DihedralType( dih->A1() + atomOffset, dih->A2() + atomOffset,
                                 dih->A3() + atomOffset, dih->A4() + atomOffset,
                                 dih->Type() ), dp[dih->Idx()] );
}

/// Map type indices to LJ parameters 
class TypeArray {
  public:
    typedef std::map<int,AtomType> Tarray;
    typedef Tarray::const_iterator const_iterator;
    typedef std::pair<int,AtomType> Tpair;
    TypeArray() {}
    unsigned int size()    const { return types_.size();  }
    const_iterator begin() const { return types_.begin(); }
    const_iterator end()   const { return types_.end();   }
    AtomType const& LastType() const { return lastType_->second; }
    /// Add params for atom to array if not present
    int AddAtomType(Atom const& atom, int anum, Topology const& top) {
      if (atom.TypeIndex() < 0 || atom.Type().len() < 1) {
        mprintf("Warning: Invalid atom type information in '%s'\n", top.c_str());
        return 1;
      }
      //mprintf("DEBUG: %s %s %i %s\n", *(atom.Name()), *(atom.Type()), anum, top.c_str());
      // See if index already present.
      Tarray::iterator it = types_.find( atom.TypeIndex() );
      if (it == types_.end()) {
        AtomType type(top.GetVDWradius(anum), top.GetVDWdepth(anum), atom.TypeIndex());
        //mprintf("\tAdding [%4zu] %4i %s %12.5g %12.5g (%s)\n", types_.size(), 
        //         atom.TypeIndex(), *(atom.Type()),
        //         type.Radius(), type.Depth(), top.c_str());
        std::pair<Tarray::iterator, bool> ret = types_.insert( Tpair(atom.TypeIndex(), type) );
        lastType_ = ret.first;
      }
      return 0;
    }
    /// Add to array
    void AddType(int idx, AtomType const& typeIn) { types_.insert(Tpair(idx, typeIn)); }
  private:
    Tarray types_;
    Tarray::iterator lastType_;
};

// Topology::AppendTop()
int Topology::AppendTop(Topology const& NewTop) {
  int atomOffset = (int)atoms_.size();
  int resOffset = (int)residues_.size();

  // NONBONDS
  // Make sure that either both topologies have non-bond parameters or neither
  // do. If the current topology is empty just check incoming topology.
  bool doNonBond = true;
  if (atoms_.empty()) {
    doNonBond = NewTop.Nonbond().HasNonbond();
  } else {
    doNonBond = (Nonbond().HasNonbond() && NewTop.Nonbond().HasNonbond());
    if (!doNonBond && (Nonbond().HasNonbond() != NewTop.Nonbond().HasNonbond())) {
      if (Nonbond().HasNonbond())
        mprintf("Warning: Topology '%s' does not have non-bond parameters.\n", NewTop.c_str());
      else
        mprintf("Warning: Topology '%s' does not have non-bond parameters.\n", c_str());
    }
  }
  // Create an array of existing nonbond types so we can compare to incoming.
  TypeArray ExistingTypes;
  if (doNonBond) {
    for (atom_iterator atom = atoms_.begin(); atom != atoms_.end(); ++atom) {
      if (ExistingTypes.AddAtomType(*atom, atom-atoms_.begin(), *this)) {
        doNonBond = false;
        break;
      }
    }
    // DEBUG
    if (debug_ > 0) {
      mprintf("DEBUG: %zu existing atom type indices:\n", ExistingTypes.size());
      for (TypeArray::const_iterator ti = ExistingTypes.begin();
                                     ti != ExistingTypes.end(); ++ti)
        mprintf("\t%8i %12.5g %12.5g\n", ti->first, ti->second.Radius(), ti->second.Depth());
    }
  }
  int NexistingAtomTypes = (int)ExistingTypes.size();
  int currentTypeIdx = NexistingAtomTypes;

  // ATOMS
  typedef std::map<int,int> TypeMap;
  TypeMap type_newToExisting; ///< Track what existing atom type new types correspond to.
  TypeArray NewTypes;
  for (atom_iterator atom = NewTop.begin(); atom != NewTop.end(); ++atom)
  {
    if (debug_ > 1)
      mprintf("DBG: %6u %s %s %4i\n", atom-NewTop.begin(), 
              *(atom->Name()), *(atom->Type()), atom->TypeIndex());
    Atom CurrentAtom = *atom;
    Residue const& res = NewTop.Res( CurrentAtom.ResNum() );
    // Bonds need to be cleared and re-added.
    CurrentAtom.ClearBonds();
    // NONBONDS
    if (doNonBond) {
      if (NewTypes.AddAtomType(*atom, atom-NewTop.begin(), NewTop)) {
        doNonBond = false;
      } else {
        // Update type index
        int newTypeIdx = -1;
        // Check if this atom type has already been added to existing types.
        TypeMap::iterator it = type_newToExisting.find( atom->TypeIndex() );
        if (it == type_newToExisting.end()) {
          // Type has not yet been mapped. See if it is an existing type.
          for (TypeArray::const_iterator et = ExistingTypes.begin();
                                         et != ExistingTypes.end(); ++et)
            if (et->second == NewTypes.LastType()) {
              newTypeIdx = et->first;
              type_newToExisting.insert(std::pair<int,int>(atom->TypeIndex(), newTypeIdx));
              if (debug_ > 1)
                mprintf("\tType (%i) matches existing type (%i)\n",
                        atom->TypeIndex(), newTypeIdx);
              break;
            }
          if (newTypeIdx == -1) {
            // Type not found among existing types. Need new type index.
            if (debug_ > 1)
              mprintf("\tNeed new type. Converting %i to %i\n",
                      atom->TypeIndex(), currentTypeIdx);
            newTypeIdx = currentTypeIdx;
            type_newToExisting.insert(std::pair<int,int>(atom->TypeIndex(), newTypeIdx));
            ExistingTypes.AddType(newTypeIdx, NewTypes.LastType());
            currentTypeIdx++;
          }
        } else {
          // Type has been mapped.
          newTypeIdx = it->second;
          if (debug_ > 1)
            mprintf("\tType %i already present as %i.\n", atom->TypeIndex(),newTypeIdx);
        }
        // Update type index
        CurrentAtom.SetTypeIndex( newTypeIdx );
      }
    }
    AddTopAtom( CurrentAtom, Residue(res.Name(), CurrentAtom.ResNum() + resOffset,
                                     res.Icode(), res.ChainID()) );
  } // END loop over incoming atoms
  // NONBONDS
  if (!doNonBond) {
    if (Nonbond().HasNonbond()) {
      mprintf("Warning: Removing non-bond parameters\n");
      nonbond_.Clear();
    }
  } else {
    // DEBUG
    if (debug_ > 0) {
      mprintf("DEBUG: %zu new atom type indices:\n", NewTypes.size());
      for (TypeArray::const_iterator ti = NewTypes.begin();
                                     ti != NewTypes.end(); ++ti)
        mprintf("\t%8i %12.5g %12.5g\n", ti->first, ti->second.Radius(), ti->second.Depth());
      mprintf("DEBUG: New to existing mapping:\n");
      for (TypeMap::const_iterator it = type_newToExisting.begin();
                                   it != type_newToExisting.end(); ++it)
        mprintf("\t%6i to %6i\n", it->first, it->second);
      mprintf("DEBUG: Atom Types:\n");
      for (TypeArray::const_iterator it = ExistingTypes.begin();
                                     it != ExistingTypes.end(); ++it)
      {
        mprintf("\t%8i %12.5g %12.5g (%8i)", it->first, it->second.Radius(),
                it->second.Depth(), it->second.OriginalIdx());
        if (it->first >= NexistingAtomTypes)
          mprintf(" (NEW)\n");
        else
          mprintf(" (EXISTING)\n");
      }
    }
    // Set up new nonbond array.
    NonbondParmType newNB;
    newNB.SetupLJforNtypes( ExistingTypes.size() );
    // Go through all interactions.
    for (TypeArray::const_iterator t1 = ExistingTypes.begin();
                                   t1 != ExistingTypes.end(); ++t1)
    {
      for (TypeArray::const_iterator t2 = ExistingTypes.begin();
                                     t2 != ExistingTypes.end(); ++t2)
      {
        if (t1->first >= t2->first)
        {
          bool OneIsNew = (t1->first >= NexistingAtomTypes);
          bool TwoIsNew = (t2->first >= NexistingAtomTypes);
          NonbondType LJ = LJ_EMPTY;
          if (OneIsNew && TwoIsNew) {
            // Both types from new topology. Look there for LJ params using
            // the original atom indices.
            int nbidx = NewTop.Nonbond().GetLJindex( t1->second.OriginalIdx(),
                                                     t2->second.OriginalIdx() );
            if (nbidx >= 0)
              LJ = NewTop.Nonbond().NBarray( nbidx );
          } else if (!OneIsNew && !TwoIsNew) {
            // Both types from existing topology. Look here for LJ params.
            int nbidx = nonbond_.GetLJindex( t1->first, t2->first );
            if (nbidx >= 0)
              LJ = nonbond_.NBarray( nbidx );
          } else {
            // Mix new and existing.
            LJ = t1->second.Combine_LB( t2->second );
          }
          if (debug_ > 1)
            mprintf("DEBUG: Adding LJ term for %i %i A=%g B=%g\n", t1->first, t2->first,
                    LJ.A(), LJ.B()); 
          newNB.AddLJterm(t1->first, t2->first, LJ);
        }
      }
    }
    nonbond_ = newNB;
  }
  // EXTRA ATOM INFO
  for (Topology::extra_iterator extra = NewTop.extraBegin();
                                extra != NewTop.extraEnd(); ++extra)
    AddExtraAtomInfo( *extra );
  // BONDS
  AddBondArray(NewTop.Bonds(),  NewTop.BondParm(), atomOffset);
  AddBondArray(NewTop.BondsH(), NewTop.BondParm(), atomOffset);
  // ANGLES
  AddAngleArray(NewTop.Angles(),  NewTop.AngleParm(), atomOffset);
  AddAngleArray(NewTop.AnglesH(), NewTop.AngleParm(), atomOffset);
  // DIHEDRALS
  AddDihArray(NewTop.Dihedrals(),  NewTop.DihedralParm(), atomOffset);
  AddDihArray(NewTop.DihedralsH(), NewTop.DihedralParm(), atomOffset);

  // Re-set up this topology
  // TODO: Could get expensive for multiple appends.
  return CommonSetup();
}
