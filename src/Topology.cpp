#include <algorithm> // find
#include <map> // For atom types in append
#include <stack> // For large system molecule search
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
#include "Constants.h" // RADDEG, SMALL
#include "AtomType.h"
#include "AtomMask.h"
#include "CharMask.h"
#include "UpdateParameters.h"

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

/** Reset all PDB-related info.
  * NOTE: This routine is used by AmbPDB.
  */
void Topology::ResetPDBinfo() {
  atom_altloc_.clear();
  occupancy_.clear();
  bfactor_.clear();
  pdbSerialNum_.clear();
  int resnum = 1;
  for (std::vector<Residue>::iterator res = residues_.begin();
                                      res != residues_.end(); ++res, ++resnum)
  {
    res->SetOriginalNum( resnum );
    res->SetIcode(' ');
    res->SetChainID(' ');
  }
}

/** Used to set box info from currently associated trajectory. */
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
    if ( boxIn.Param(Box::X) < Constants::SMALL || 
         boxIn.Param(Box::Y) < Constants::SMALL || 
         boxIn.Param(Box::Z) < Constants::SMALL )
    {
      // Incoming box has no lengths - disable parm box. TODO is this check necessary/desirable?
      mprintf("Warning: Box information present in trajectory but lengths are zero.\n"
              "Warning: DISABLING BOX in topology '%s'!\n", c_str());
      parmBox_.SetNoBox();
    } else {
      // Incoming box is valid. Indicate if current box type differs from
      // incoming box type.
      if (parmBox_.CellShape() != boxIn.CellShape()) {
        mprintf("Warning: Trajectory box type is '%s' but topology box type is '%s'.\n"
                "Warning: Setting topology box information from trajectory.\n",
                boxIn.CellShapeName(), parmBox_.CellShapeName());
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
  if (molecules_.size() > 0) {
    // Topology has molecule information
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
  } else {
    // No molecule information
    mprintf("Warning: No molecule information. Determining solvent residues based on naming.\n");
    for (int res = 0; res != Nres(); res++) {
      Residue const& currentRes = Res(res);
      if (!currentRes.NameIsSolvent()) {
        // Not a solvent name.
        if (currentRes.NumAtoms() > 1 ||
            Atoms()[currentRes.FirstAtom()].BondIdxArray().size() > 0)
        {
          // If this residue has > 1 atom, or is 1 atom but has bonds, assume solute.
          solute_res.AddToRange( res );
        }
      }
    } // END loop over residues
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

/** Given an atom number, return a string containing the corresponding
  * residue name and atom name with format:
  * "<resname>@<atom name>"
  * Truncate the residue and atom names so there are no blanks.
  */
std::string Topology::TruncResNameAtomName(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string("");
  int res = atoms_[atom].ResNum();
  return residues_[res].Name().Truncated() + "@" + atoms_[atom].Name().Truncated();
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

/** Given an atom number, return a string containing the corresponding
  * residue name and number, and atom name and number, all separated
  * by spaces:
  * "<resname> <resnum> <atom name> <atom num>
  */
std::string Topology::ResNameNumAtomNameNum(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string("");
  int res = atoms_[atom].ResNum();
  return residues_[res].Name().Truncated() + " " + integerToString(res+1) + " " +
         atoms_[atom].Name().Truncated() + " " + integerToString(atom+1);
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

/** Given an atom number, return a string containing atom name and
  * number with format:
  * "<atomname>_<atomnum>"
  */
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
std::string Topology::TruncResNameNum(int res) const {
  if (res < 0 || res >= (int)residues_.size()) return std::string("");
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
  mprintf("\t\tBox: %s\n", parmBox_.CellShapeName());
  if (NsolventMolecules_>0) {
    mprintf("\t\t%i solvent molecules.\n", NsolventMolecules_);
  }
  if (!radius_set_.empty())
    mprintf("\t\tGB radii set: %s\n", radius_set_.c_str());
  if (nonbond_.HasNonbond()) {
    mprintf("\t\tNon-bonded parameters are present.\n");
    if (nonbond_.Has_C_Coeff())
      mprintf("\t\t\tLJ 12-6-4 C coefficients are present.\n");
  }
  if (chamber_.HasChamber()) {
    mprintf("\t\tCHAMBER: %zu Urey-Bradley terms, %zu Impropers\n",
            chamber_.UB().size(), chamber_.Impropers().size());
    if (HasCmap())
      mprintf("\t\t         %zu CMAP grids, %zu CMAP terms.\n", 
              CmapGrid().size(), Cmap().size());
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
          residues_.size(), parmBox_.CellShapeName(), molecules_.size());
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
  if ( residues_.empty() || residues_.back() != resIn )
  {
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
/*
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
}*/

/** Check the given size. If not # atoms, return true (error). */
bool Topology::CheckExtraSize(size_t sizeIn, const char* desc)
const
{
  if (sizeIn > 0 && sizeIn != atoms_.size()) {
    mprinterr("Error: Size of the %s array (%zu) is not # atoms (%zu)\n", desc, sizeIn, atoms_.size());
    return true;
  }
  return false;
}

// Topology::CommonSetup()
/** Set up common to all topologies.
  * \param molsearch If true, determine molecules based on bond info.
  * \param renumberResidues If true, renumber residues if any residue is part of more than 1 molecule
  *        e.g. when alternate locations are present.
  */
int Topology::CommonSetup(bool molsearch, bool renumberResidues)
{
  // Check the size of any "extra" arrays
  if (CheckExtraSize(tree_.size(), "Amber tree")) return 1;
  if (CheckExtraSize(ijoin_.size(), "Amber join")) return 1;
  if (CheckExtraSize(irotat_.size(), "Amber rotate")) return 1;
  if (CheckExtraSize(atom_altloc_.size(), "PDB alt. loc.")) return 1;
  if (CheckExtraSize(occupancy_.size(), "PDB occupancy")) return 1;
  if (CheckExtraSize(bfactor_.size(), "PDB Bfactor")) return 1;
  if (CheckExtraSize(pdbSerialNum_.size(), "PDB serial #")) return 1;
  // TODO: Make bond parm assignment / molecule search optional?
  // Assign default lengths if necessary (for e.g. CheckStructure)
  if (bondparm_.empty())
    AssignBondParameters();
  if (molsearch) {
    // Determine molecule info from bonds
    if (DetermineMolecules())
      mprinterr("Error: Could not determine molecule information for %s.\n", c_str());
  }
  // DEBUG : Current residue info
  if (debug_ > 1) {
    mprintf("DEBUG: Current residue info (%zu).\n", residues_.size());
    for (std::vector<Residue>::const_iterator res = residues_.begin(); res != residues_.end(); ++res)
    {
      mprintf("DEBUG:\t\t%8li %6s orig=%8i atoms %8i to %8i\n", res-residues_.begin(),
              *(res->Name()), res->OriginalResNum(), res->FirstAtom(), res->LastAtom());
    }
  }
  // Check if any molecules share residue numbers. If so and if specified,
  // base residue information on molecules.
  if (renumberResidues && !molecules_.empty() && molecules_.size() > 1) {
    bool mols_share_residues = (molecules_.size() > residues_.size());
    if (!mols_share_residues) {
      // More in-depth check
      for (std::vector<Molecule>::const_iterator mol = molecules_.begin() + 1;
                                                 mol != molecules_.end(); ++mol)
      {
        int m0_resnum = atoms_[(mol-1)->MolUnit().Front()].ResNum();
        int m1_resnum = atoms_[    mol->MolUnit().Front()].ResNum();
        if (m0_resnum == m1_resnum) {
          mols_share_residues = true;
          long int molnum = mol - molecules_.begin();
          mprintf("Warning: 2 or more molecules (%li and %li) share residue numbers (%i).\n",
                  molnum, molnum+1, m0_resnum+1);
          break;
        }
      }
    }
    if (mols_share_residues) {
      mprintf("Warning:   This usually happens when alternate locations for atoms are present.\n"
              "Warning:   Basing residue information on molecules.\n");
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
  } // END renumber residues based on molecules

  // Set up solvent information
  if (SetSolventInfo())
    mprinterr("Error: Could not determine solvent information for %s.\n", c_str());

  // Determine # of extra points.
  DetermineNumExtraPoints();

  return 0;
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
    if (mol->MolUnit().nSegments() == 1 && mol->NumAtoms() == 3) {
      int nH = 0;
      int nO = 0;
      for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                seg != mol->MolUnit().segEnd(); ++seg)
      {
        for (int atnum = seg->Begin(); atnum != seg->End(); atnum++)
        {
          if (atoms_[atnum].Element() == Atom::HYDROGEN) nH++;
          if (atoms_[atnum].Element() == Atom::OXYGEN)   nO++;
        }
      }
      if (nO == 1 && nH == 2) res_name = "HOH";
    } else
      res_name = default_res_name;
    residues_.push_back( Residue(res_name, resnum+1, ' ', ' ') );
    residues_.back().SetFirstAtom( mol->MolUnit().Front() );
    residues_.back().SetLastAtom( mol->MolUnit().Back() );
    // Update atom residue numbers
    for (int atnum = residues_.back().FirstAtom(); 
             atnum != residues_.back().LastAtom(); ++atnum)
      atoms_[atnum].SetResNum( resnum );
  }
  return 0;
}

// Topology::Resize()
/** Clear all arrays; allocate atoms, residues, tree, ijoin, irotat, and
  * bond/angle/dihedral parameter arrays according to input pointers.
  * Intended for use when reading Amber Topology file, specifically the
  * POINTERS section.
  */
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
  cmap_.clear();
  cmapGrid_.clear();
  nonbond_.Clear();
  cap_.Clear();
  lesparm_.Clear();
  chamber_.Clear();
  tree_.clear();
  ijoin_.clear();
  irotat_.clear();
  atom_altloc_.clear();
  occupancy_.clear();
  bfactor_.clear();
  pdbSerialNum_.clear();
  parmBox_.SetNoBox();
  refCoords_ = Frame();
  ipol_ = 0;
  NsolventMolecules_ = 0;
  n_extra_pts_ = 0;
  n_atom_types_ = 0;

  atoms_.resize( pIn.natom_ );
  residues_.resize( pIn.nres_ );
  tree_.resize( pIn.natom_ );
  ijoin_.resize( pIn.natom_, 0 );
  irotat_.resize( pIn.natom_, 0 );
  bondparm_.resize( pIn.nBndParm_ );
  angleparm_.resize( pIn.nAngParm_ );
  dihedralparm_.resize( pIn.nDihParm_ );
}

/** \return Rmin for given atom. */
double Topology::GetVDWradius(int a1) const {
  //TODO: return zero when no params?
  return GetLJparam(a1, a1).Radius();
}

/** \return sigma for given atom. */
double Topology::GetVDWsigma(int a1) const {
  //TODO: return zero when no params?
  NonbondType const& LJ = GetLJparam(a1, a1);
  if (LJ.B() > 0.0)
    return ( 0.5 * pow(LJ.A() / LJ.B(), (1.0/6.0)) );
  else
    return 0.0;
}

/** \return epsilon for given atom. */
double Topology::GetVDWdepth(int a1) const {
  return GetLJparam(a1, a1).Depth();
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
  mprintf("Warning: Determining bond length parameters from element types for '%s'.\n", c_str());
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

/** Remove a bond between atom 1 and atom2, update the atoms array.
  * Does not modify bond parameters.
  * \return 0 if a bond was successfully removed, -1 if no bond exists, and 1 if an error occurs.
  */
int Topology::RemoveBond(int atom1, int atom2)
{
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), atom1, "bond")) return 1;
  if (WarnOutOfRange(atoms_.size(), atom2, "bond")) return 1;
  // Ensure the bond exists.
  bool exists = false;
  for (Atom::bond_iterator ba = atoms_[atom1].bondbegin();
                           ba != atoms_[atom1].bondend(); ++ba)
    if ( *ba == atom2 ) {
      exists = true;
      break;
    }
  if (!exists) {
    mprintf("Warning: No bond exists between atoms %i and %i\n", atom1+1, atom2+1);
    return -1;
  }
  bool a1H = (atoms_[atom1].Element() == Atom::HYDROGEN);
  bool a2H = (atoms_[atom2].Element() == Atom::HYDROGEN);
  BondArray* tgtArray;
  if (a1H || a2H)
    tgtArray = &bondsh_;
  else
    tgtArray = &bonds_;
  // Search the array.
  BondArray::iterator bnd = tgtArray->begin();
  for (; bnd != tgtArray->end(); ++bnd) {
    if (atom1 == bnd->A1()) {
      if (atom2 == bnd->A2()) break;
    }
    if (atom2 == bnd->A1()) {
      if (atom1 == bnd->A2()) break;
    }
  }
  // Sanity check
  if (bnd == tgtArray->end()) {
    mprinterr("Internal Error: Bond %i %i not found in internal bond array.\n", atom1+1, atom2+1);
    return 1;
  }
  tgtArray->erase( bnd );
  atoms_[atom1].RemoveBondToIdx( atom2 );
  atoms_[atom2].RemoveBondToIdx( atom1 );
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

// -----------------------------------------------
// Topology::AddTorsionParm()
/** Check if given dihedral parm exists in given dihedral parm array. Add if not.
  * \return Index in dihedral parm array.
  */
int Topology::AddTorsionParm(DihedralParmArray& dparray, DihedralParmType const& DPin)
{
  // See if the DihedralParm exists.
  int pidx = -1;
  for (DihedralParmArray::const_iterator dp = dparray.begin();
                                         dp != dparray.end(); ++dp)
    if ( fabs(DPin.Pk()    - dp->Pk()   ) < Constants::SMALL &&
         fabs(DPin.Pn()    - dp->Pn()   ) < Constants::SMALL &&
         fabs(DPin.Phase() - dp->Phase()) < Constants::SMALL &&
         fabs(DPin.SCEE()  - dp->SCEE() ) < Constants::SMALL &&
         fabs(DPin.SCNB()  - dp->SCNB() ) < Constants::SMALL )
    {
      pidx = (int)(dp - dparray.begin());
      break;
    }
  if (pidx == -1) {
    pidx = (int)dparray.size();
    dparray.push_back( DPin );
  }
  return pidx;
}

/** \return true if any atoms in the dihedral are out of range. */
bool Topology::CheckTorsionRange(DihedralType const& dihIn, const char* typestr) const
{
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), dihIn.A1(), typestr)) return true;
  if (WarnOutOfRange(atoms_.size(), dihIn.A2(), typestr)) return true;
  if (WarnOutOfRange(atoms_.size(), dihIn.A3(), typestr)) return true;
  if (WarnOutOfRange(atoms_.size(), dihIn.A4(), typestr)) return true;
  return false;
}

/** \return Dihedral with parm index set. */
DihedralType Topology::SetTorsionParmIndex(DihedralType const& dihIn,
                                           DihedralParmArray const& dparray,
                                           int pidxIn, const char* typestr)
{
  // Check if parm index is out of range;
  int pidx;
  if (pidxIn < (int)dparray.size())
    pidx = pidxIn;
  else {
    mprintf("Warning: No %s parameters for index %i\n", typestr, pidxIn);
    pidx = -1;
  }
  DihedralType dih = dihIn;
  dih.SetIdx( pidx );
  return dih;
}

/** Add given dihedral with given dihedral parm to dihedral array. */
void Topology::AddDihedral(DihedralType const& dih, DihedralParmType const& DPin)
{
  int pidx = AddTorsionParm(dihedralparm_, DPin);
  if (CheckTorsionRange(dih, "dihedral")) return;
  AddDihedral(dih, pidx);
}

/** Add given dihedral with given dihedral parm index. */
void Topology::AddDihedral(DihedralType const& dihIn, int pidxIn) {
  // FIXME: Check duplicate
  if (CheckTorsionRange(dihIn, "dihedral")) return;
  DihedralType dih = SetTorsionParmIndex(dihIn, dihedralparm_, pidxIn, "dihedral");
  // Update dihedral arrays
  if (atoms_[dih.A1()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A2()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A3()].Element() == Atom::HYDROGEN ||
      atoms_[dih.A4()].Element() == Atom::HYDROGEN)
    dihedralsh_.push_back( dih );
  else
    dihedrals_.push_back( dih );
}

/** Add given dihedral to either dihedral with or without H array. */
void Topology::AddDihedral(DihedralType const& dihIn, bool isH) {
  if (isH)
    dihedralsh_.push_back( dihIn );
  else
    dihedrals_.push_back( dihIn );
}

/** Add given Charmm improper with given improper parm to Charmm improper array. */
void Topology::AddCharmmImproper(DihedralType const& imp, DihedralParmType const& IPin)
{
  int pidx = AddTorsionParm(chamber_.SetImproperParm(), IPin);
  if (CheckTorsionRange(imp, "CHARMM improper")) return;
  AddCharmmImproper(imp, pidx);
}

/** Add given Charmm improper with given improper parm index. */
void Topology::AddCharmmImproper(DihedralType const& impIn, int pidxIn)
{
  if (CheckTorsionRange(impIn, "CHARMM improper")) return;
  DihedralType imp = SetTorsionParmIndex(impIn, chamber_.ImproperParm(), pidxIn, "CHARMM improper");
  // Update Charmm improper array.
  chamber_.AddImproperTerm( imp );
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

/** \return Number of residues in specified molecule. */
int Topology::NresInMol(int idx) const {
  int nres = 0;
  for (Unit::const_iterator seg = molecules_[idx].MolUnit().segBegin();
                            seg != molecules_[idx].MolUnit().segEnd(); ++seg)
    nres += atoms_[seg->End()-1].ResNum() - atoms_[seg->Begin()].ResNum() + 1;
  return nres;
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
  if (numberOfMolecules < 1) {
    mprinterr("Internal Error: Could not determine molecules.\n");
    return 1;
  }
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
      mprintf("\t\tAtom %li assigned to molecule %i\n", atom - atoms_.begin(), atom->MolNum());
  }

  // Update molecule information
  molecules_.resize( numberOfMolecules );
  for (int atomIdx = 0; atomIdx < (int)atoms_.size(); atomIdx++)
  {
    Atom const& atom = atoms_[atomIdx];
    molecules_[atom.MolNum()].ModifyUnit().AddIndex( atomIdx );
  }
  if (debug_ > 0) mprintf("DEBUG: Molecule segment information:\n");
  std::vector< std::vector<Molecule>::const_iterator > nonContiguousMols;
  for (std::vector<Molecule>::const_iterator mol = molecules_.begin(); mol != molecules_.end(); ++mol)
  {
    if (mol->MolUnit().nSegments() > 1)
      nonContiguousMols.push_back( mol );
    if (debug_ > 0) {
      mprintf("DEBUG:\t%8li %8u segments:", mol - molecules_.begin() + 1, mol->MolUnit().nSegments());
      for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                         seg != mol->MolUnit().segEnd(); ++seg)
        mprintf(" %i-%i (%i) ", seg->Begin()+1, seg->End(), seg->Size());
      mprintf("\n");
    }
  }
  if (!nonContiguousMols.empty()) {
    mprintf("Warning: %zu molecules have non-contiguous segments of atoms.\n", nonContiguousMols.size());
    for (std::vector< std::vector<Molecule>::const_iterator >::const_iterator it = nonContiguousMols.begin();
                                                                              it != nonContiguousMols.end(); ++it)
    {
      mprintf("\t%8li %8u segments:", *it - molecules_.begin() + 1, (*it)->MolUnit().nSegments());
      for (Unit::const_iterator seg = (*it)->MolUnit().segBegin();
                                seg != (*it)->MolUnit().segEnd(); ++seg)
        mprintf(" %i-%i (%i) ", seg->Begin()+1, seg->End(), seg->Size());
      mprintf("\n");
    }
    mprintf("Warning: The 'fixatomorder' command can be used to reorder the topology and any\n"
            "Warning:  associated coordinates.\n");
  } 
/*
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
      mprinterr("Error: Atom %li was assigned a lower molecule # (%i) than previous atom (%i).\n"
                "Error:   This can happen if bond information is incorrect or missing, or if the\n"
                "Error:   atom numbering in molecules is not sequential. Try one of the\n"
                "Error:   following:\n"
                "Error: - If this is a PDB file, try using the 'noconect' keyword.\n"
                "Error: - If this topology did not have bond info, try increasing the bond\n"
                "Error:   search cutoff above 0.2 Ang. ('bondsearch <cutoff>').\n"
                "Error: - Use the 'fixatomorder' command to reorder the topology and any\n"
                "Error:   associated coordinates.\n"
                "Error: - Use the 'setMolecules' command in parmed to reorder only the\n"
                "Error:   topology.\n", atom - atoms_.begin() + 1,
                atom->MolNum()+1, lastMol+1);
      ClearMolecules();
      return 1;
    }
    ++atomNum;
  }
  molecule->SetLast( atoms_.size() );
*/
  return 0;
}

// -----------------------------------------------------------------------------
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
  mask.MaskInfo();
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
    if ( mask.AtomsInCharMask( mol->MolUnit() ) ) {
      mol->SetSolvent();
      ++NsolventMolecules_;
      numSolvAtoms += mol->NumAtoms();
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
    int firstRes = atoms_[ mol->MolUnit().Front() ].ResNum();
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

//  Topology::ResnumsSelectedBy()
std::vector<int> Topology::ResnumsSelectedBy(AtomMask const& mask) const {
  std::vector<int> resnums;
  int res = -1;
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
    if (atoms_[*at].ResNum() > res) {
      res = atoms_[*at].ResNum();
      resnums.push_back( res );
    }
  return resnums;
}

// Topology::MolnumsSelectedBy()
std::vector<int> Topology::MolnumsSelectedBy(AtomMask const& mask) const {
  std::set<int> molnums;
  if (molecules_.empty()) {
    mprintf("Warning: Topology has no molecule information.\n");
  } else {
    int mol = -1;
    for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
      if (atoms_[*at].MolNum() != mol) {
        mol = atoms_[*at].MolNum();
        molnums.insert( mol );
      }
  }
  std::vector<int> tmp;
  tmp.reserve( molnums.size() );
  for (std::set<int>::const_iterator it = molnums.begin(); it != molnums.end(); ++it)
    tmp.push_back( *it );
  return tmp;
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

/** This template can be used when doing ModifyByMap() on a generic std::vector array
  * of type T. 
  */
template <class T> class TopVecStrip {
  public:
    /// CONSTRUCTOR
    TopVecStrip() {}
    /// Strip current array according to given map, output to given array of same type
    void Strip(std::vector<T>& newArray, std::vector<T> const& oldArray, std::vector<int> const& MapIn)
    {
      if (!oldArray.empty()) {
        for (std::vector<int>::const_iterator old_it = MapIn.begin(); old_it != MapIn.end(); ++old_it)
          if (*old_it >= 0)
            newArray.push_back( oldArray[*old_it] );
      }
    }
};

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
      newParm->residues_.push_back( cr );
      newParm->residues_.back().SetFirstAtom( newatom );
      newParm->residues_.back().SetTerminal( cr.IsTerminal() );
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
      if ( isSolvent[ mol->MolUnit().Front() ] ) {
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
    if (HasCmap()) {
      // NOTE that atom indexing is updated but cmap indexing is not. So if
      // any CMAP terms remain all CMAP entries remain.
      for (CmapArray::const_iterator cmap = Cmap().begin();
                                     cmap != Cmap().end(); ++cmap)
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
                  newParm->AddCmapTerm( CmapType(newA1,newA2,newA3,
                                                          newA4,newA5,cmap->Idx()) );
              }
            }
          }
        }
      }
      // Only add CMAP grids if there are CMAP terms left.
      if (!newParm->Cmap().empty()) {
        for (CmapGridArray::const_iterator g = CmapGrid().begin();
                                           g != CmapGrid().end(); ++g)
          newParm->AddCmapGrid( *g );
      }
    }
  }
  // Amber extra info.
  TopVecStrip<NameType> stripNameType;
  stripNameType.Strip(newParm->tree_, tree_, MapIn);
  TopVecStrip<int> stripInt;
  stripInt.Strip(newParm->ijoin_, ijoin_, MapIn);
  stripInt.Strip(newParm->irotat_, irotat_, MapIn);
  stripInt.Strip(newParm->pdbSerialNum_, pdbSerialNum_, MapIn);
  TopVecStrip<char> stripChar;
  stripChar.Strip(newParm->atom_altloc_, atom_altloc_, MapIn);
  TopVecStrip<float> stripFloat;
  stripFloat.Strip(newParm->occupancy_, occupancy_, MapIn);
  stripFloat.Strip(newParm->bfactor_, bfactor_, MapIn);
  
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
    if (oldidx > -1) {
      int newidx = parmMap[oldidx];
      if (newidx == -1) { // This needs to be added to new parameter array.
        newidx = (int)newBondParm.size();
        parmMap[oldidx] = newidx;
        newBondParm.push_back( oldParm[oldidx] );
      }
      //mprintf("DEBUG: Old bond parm index=%i, new bond parm index=%i\n", oldidx, newidx);
      bnd->SetIdx( newidx );
    }
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
    if (oldidx > -1) {
      int newidx = parmMap[oldidx];
      if (newidx == -1) { // This needs to be added to new parameter array.
        newidx = (int)newAngleParm.size();
        parmMap[oldidx] = newidx;
        newAngleParm.push_back( angleparm_[oldidx] );
      }
      //mprintf("DEBUG: Old angle parm index=%i, new angle parm index=%i\n", oldidx, newidx);
      ang->SetIdx( newidx );
    }
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
    if (oldidx > -1) {
      int newidx = parmMap[oldidx];
      if (newidx == -1) { // This needs to be added to new parameter array.
        newidx = (int)newDihedralParm.size();
        parmMap[oldidx] = newidx;
        newDihedralParm.push_back( oldParm[oldidx] );
      }
      //mprintf("DEBUG: Old dihedral parm index=%i, new dihedral parm index=%i\n", oldidx, newidx);
      dih->SetIdx( newidx );
    }
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
    size_t size()    const { return types_.size();  }
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

/** This template can be used when doing Append() on a generic std::vector array
  * of type T. The array will be appended to a given array of the same type.
  * If one is empty and the other is not, values will be filled in if necessary.
  */
template <class T> class TopVecAppend {
  public:
    /// CONSTRUCTOR
    TopVecAppend() {}
    /// Append current array to given array of same type
    void Append(std::vector<T>& arrayOut, std::vector<T> const& arrayToAdd, unsigned int expectedSize)
    {
      if (arrayToAdd.empty() && arrayOut.empty()) {
        // Both arrays are empty. Nothing to do.
        return;
      } else if (arrayToAdd.empty()) {
        // The current array is empty but the given array is not. Fill in 
        // array to append with blank values.
        for (unsigned int idx = 0; idx != expectedSize; idx++)
          arrayOut.push_back( T() );
      } else {
        // Append current array to array to given array. TODO use std::copy?
        for (typename std::vector<T>::const_iterator it = arrayToAdd.begin(); it != arrayToAdd.end(); ++it)
          arrayOut.push_back( *it );
      }
    }
};

// Topology::AppendTop()
int Topology::AppendTop(Topology const& NewTop) {
  int atomOffset = (int)atoms_.size();
  //int resOffset = (int)residues_.size();

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
        mprintf("\t%8i %12.5g %12.5g\n", ti->first, ti->second.LJ().Radius(), ti->second.LJ().Depth());
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
      mprintf("DBG: %6li %s %s %4i\n", atom-NewTop.begin(), 
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
    //AddTopAtom( CurrentAtom, Residue(res.Name(), CurrentAtom.ResNum() + resOffset + 1,
    AddTopAtom( CurrentAtom, Residue(res.Name(), res.OriginalResNum(),
                                     res.Icode(), res.ChainId()) );
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
        mprintf("\t%8i %12.5g %12.5g\n", ti->first, ti->second.LJ().Radius(), ti->second.LJ().Depth());
      mprintf("DEBUG: New to existing mapping:\n");
      for (TypeMap::const_iterator it = type_newToExisting.begin();
                                   it != type_newToExisting.end(); ++it)
        mprintf("\t%6i to %6i\n", it->first, it->second);
      mprintf("DEBUG: Atom Types:\n");
      for (TypeArray::const_iterator it = ExistingTypes.begin();
                                     it != ExistingTypes.end(); ++it)
      {
        mprintf("\t%8i %12.5g %12.5g (%8i)", it->first, it->second.LJ().Radius(),
                it->second.LJ().Depth(), it->second.OriginalIdx());
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
            LJ = t1->second.LJ().Combine_LB( t2->second.LJ() );
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
  TopVecAppend<NameType> appendNameType;
  appendNameType.Append( tree_, NewTop.tree_, NewTop.Natom() );
  TopVecAppend<int> appendInt;
  appendInt.Append( ijoin_, NewTop.ijoin_, NewTop.Natom() );
  appendInt.Append( irotat_, NewTop.irotat_, NewTop.Natom() );
  appendInt.Append( pdbSerialNum_, NewTop.pdbSerialNum_, NewTop.Natom() );
  TopVecAppend<char> appendChar;
  appendChar.Append( atom_altloc_, NewTop.atom_altloc_, NewTop.Natom() );
  TopVecAppend<float> appendFloat;
  appendFloat.Append( occupancy_, NewTop.occupancy_, NewTop.Natom() );
  appendFloat.Append( bfactor_, NewTop.bfactor_, NewTop.Natom() );

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

// -----------------------------------------------------------------------------
static void paramOverwriteWarning(const char* type) {
  mprintf("Warning: An existing %s parameter would have been overwritten. This\n"
          "Warning:  usually means that the atom type information in the Topology is\n"
          "Warning:  incomplete. This can happen for example with Chamber topologies\n"
          "Warning:  if the original atom type names were > 4 characters.\n", type);
  mprintf("Warning: The %s parameters in this topology may now be incorrect.\n", type);
}

// GetBondParams()
static inline void GetBondParams(ParmHolder<BondParmType>& BP, std::vector<Atom> const& atoms, BondArray const& bonds, BondParmArray const& bpa) {
  for (BondArray::const_iterator b = bonds.begin(); b != bonds.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(2);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      ParameterHolders::RetType ret = BP.AddParm( types, bpa[b->Idx()], false );
      if (ret == ParameterHolders::ERR)
        paramOverwriteWarning("bond");
    }
  }
}

// GetAngleParams()
static inline void GetAngleParams(ParmHolder<AngleParmType>& AP, std::vector<Atom> const& atoms, AngleArray const& angles, AngleParmArray const& apa) {
  for (AngleArray::const_iterator b = angles.begin(); b != angles.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(3);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      ParameterHolders::RetType ret = AP.AddParm( types, apa[b->Idx()], false );
      if (ret == ParameterHolders::ERR)
        paramOverwriteWarning("angle");
    }
  }
}

// GetImproperParams()
static inline void GetImproperParams(ParmHolder<DihedralParmType>& IP, std::vector<Atom> const& atoms, DihedralArray const& imp, DihedralParmArray const& ipa) {
  for (DihedralArray::const_iterator b = imp.begin(); b != imp.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(4);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      types.AddName( atoms[b->A4()].Type() );
      ParameterHolders::RetType ret = IP.AddParm( types, ipa[b->Idx()], false );
      if (ret == ParameterHolders::ERR)
        paramOverwriteWarning("improper");
    }
  }
}

// GetDihedralParams()
static inline void GetDihedralParams(DihedralParmHolder& DP, std::vector<Atom> const& atoms, DihedralArray const& dih, DihedralParmArray const& dpa) {
  for (DihedralArray::const_iterator b = dih.begin(); b != dih.end(); ++b)
  {
    if (b->Idx() != -1) {
      TypeNameHolder types(4);
      types.AddName( atoms[b->A1()].Type() );
      types.AddName( atoms[b->A2()].Type() );
      types.AddName( atoms[b->A3()].Type() );
      types.AddName( atoms[b->A4()].Type() );
      //mprintf("DEBUG: dihedral %li ( %i %i %i %i )\n", b - dih.begin() + 1, b->A1()+1, b->A2()+1, b->A3()+1, b->A4()+1);
      ParameterHolders::RetType ret = DP.AddParm( types, dpa[b->Idx()], false );
      if (ret == ParameterHolders::ERR)
        paramOverwriteWarning("dihedral");
    }
  }
}

/** Used to map type names to type indices to access nonbond parameters. */
typedef std::map<NameType, int> NameIdxMapType;

/** \param atomTypes Output array of atom types.
  * \param NB1 Current nonbond parameters.
  * \param atoms Array of atoms.
  * \param NB0 Output array of nonbond parameters.
  */
static inline void GetLJAtomTypes( ParmHolder<AtomType>& atomTypes, ParmHolder<NonbondType>& NB1, std::vector<Atom> const& atoms, NonbondParmType const& NB0, int debugIn)
{
  // TODO check for off-diagonal terms
  if (NB0.HasNonbond()) {
    // Map type names to type indices to access nonbond parameters.
    NameIdxMapType nameIdxMap;
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
    {
      // TODO check for blank type name?
      TypeNameHolder types(1);
      types.AddName( atm->Type() );
      int idx = NB0.GetLJindex( atm->TypeIndex(), atm->TypeIndex() );
      ParameterHolders::RetType ret;
      if (idx > -1) {
        NonbondType const& LJ = NB0.NBarray( idx );
        ret = atomTypes.AddParm( types, AtomType(LJ.Radius(), LJ.Depth(), atm->Mass()), true );
      } else
        ret = atomTypes.AddParm( types, AtomType(atm->Mass()), true );
      if (ret == ParameterHolders::ADDED)
        nameIdxMap.insert( std::pair<NameType, int>(atm->Type(), atm->TypeIndex()) );
    }
    // Do atom type pairs
    unsigned int nModifiedOffDiagonal = 0;
    for (ParmHolder<AtomType>::const_iterator i1 = atomTypes.begin(); i1 != atomTypes.end(); ++i1)
    {
      for (ParmHolder<AtomType>::const_iterator i2 = i1; i2 != atomTypes.end(); ++i2)
      {
        // Determine what A and B parameters would be.
        NameType const& name1 = i1->first[0];
        NameType const& name2 = i2->first[0];
        AtomType const& type1 = i1->second;
        AtomType const& type2 = i2->second;
        NonbondType lj0 = type1.LJ().Combine_LB( type2.LJ() );
        // Extract original A and B parameters.
        NameIdxMapType::const_iterator t1 = nameIdxMap.find( name1 );
        NameIdxMapType::const_iterator t2 = nameIdxMap.find( name2 );
        int idx1 = t1->second;
        int idx2 = t2->second;
        int idx = NB0.GetLJindex( idx1, idx2 );
        if (idx < 0) {
          mprinterr("Error: No off-diagonal LJ for  %s %s (%i %i)\n",
                    *name1, *name2, idx1, idx2);
          return;
        }
        NonbondType lj1 = NB0.NBarray( idx );
        // Compare them
        if (lj0 != lj1) {
          nModifiedOffDiagonal++;
          if (debugIn > 0) {
            double deltaA = fabs(lj0.A() - lj1.A());
            double deltaB = fabs(lj0.B() - lj1.B());
            mprintf("DEBUG: Potential off-diagonal LJ: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                    *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
            mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
          }
        }
        TypeNameHolder types(2);
        types.AddName( name1 );
        types.AddName( name2 );
        NB1.AddParm( types, lj1, false );
      } // END inner loop over atom types
    } // END outer loop over atom types
    if (nModifiedOffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ terms present.\n", nModifiedOffDiagonal);
  } else {
    for (std::vector<Atom>::const_iterator atm = atoms.begin(); atm != atoms.end(); ++atm)
      if (atm->Type().len() > 0)
        atomTypes.AddParm( TypeNameHolder(atm->Type()), AtomType(atm->Mass()), true );
  }
}

/** \return ParameterSet for this Topology. */
ParameterSet Topology::GetParameters() const {
  ParameterSet Params;
  // Atom LJ types
  GetLJAtomTypes( Params.AT(), Params.NB(), atoms_, nonbond_, debug_ );
  // Bond parameters.
  GetBondParams( Params.BP(), atoms_, bonds_, bondparm_ );
  GetBondParams( Params.BP(), atoms_, bondsh_, bondparm_ );
  // Angle parameters.
  GetAngleParams( Params.AP(), atoms_, angles_, angleparm_);
  GetAngleParams( Params.AP(), atoms_, anglesh_, angleparm_);
  // Dihedral parameters.
  GetDihedralParams( Params.DP(), atoms_, dihedrals_, dihedralparm_);
  GetDihedralParams( Params.DP(), atoms_, dihedralsh_, dihedralparm_);
  // CHARMM parameters
  if (chamber_.HasChamber()) {
    // UB parameters
    GetBondParams(Params.UB(), atoms_, chamber_.UB(), chamber_.UBparm());
    // Impropers
    GetImproperParams( Params.IP(), atoms_, chamber_.Impropers(), chamber_.ImproperParm() );
  }

  return Params;
}

// -----------------------------------------------------------------------------
/** Set parameters in array V for objects in array U using parameters from array T. 
  */
/*
template<typename T, typename U, typename V>
void AssignParm(std::vector<Atom> const& atoms,
                ParmHolder<T> const& newParams,  // BondParmType
                ParmHolder<int>& currentIndices,
                U& objects,                      // BondArray
                V& currentParams)                // BondParmArray
{
  for (typename U::iterator obj = objects.begin(); obj != objects.end(); ++obj) {
    TypeNameHolder types(obj->Nidx());
    for (unsigned int idx = 0; idx != obj->Nidx(); idx++)
      types.AddName( atoms[obj->Atom(idx)].Type() );
    bool found;
    // See if parameter already present.
    int idx = currentIndices.FindParam( types, found );
    if (!found) {
      // Search in new
      T param = newParams.FindParam( types, found );
      if (found) {
        // Add parameter
        idx = (int)currentParams.size();
        currentParams.push_back( param );
      } else
        idx = -1;
    }
    //if (idx == -1)
    //  mprintf("Warning: Bond parameter not found for bond %s-%s (%s-%s)\n",
    //          TruncResAtomNameNum(bnd->A1()).c_str(),
    //          TruncResAtomNameNum(bnd->A2()).c_str(),
    //          *types[0], *types[1]);
    obj->SetIdx( idx );
  }
}*/

/** Set parameters for bonds in given bond array. */
void Topology::AssignBondParm(ParmHolder<BondParmType> const& newBondParams,
                              ParmHolder<int>& currentIndices,
                              BondArray& bonds, BondParmArray& bpa, const char* desc)
{
  for (BondArray::iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    TypeNameHolder types(2);
    types.AddName( atoms_[bnd->A1()].Type() );
    types.AddName( atoms_[bnd->A2()].Type() );
    bool found;
    // See if parameter already present.
    int idx = currentIndices.FindParam( types, found );
    if (!found) {
      // Search in new
      BondParmType bp = newBondParams.FindParam( types, found );
      if (found) {
        // Add parameter
        idx = (int)bpa.size();
        bpa.push_back( bp );
        currentIndices.AddParm( types, idx, false );
      } else
        idx = -1;
    }
    if (idx == -1)
      mprintf("Warning: parameter not found for %s %s-%s (%s-%s)\n", desc,
              TruncResAtomNameNum(bnd->A1()).c_str(),
              TruncResAtomNameNum(bnd->A2()).c_str(),
              *types[0], *types[1]);
    bnd->SetIdx( idx );
  }
}

/** Replace any current bond parameters with given bond parameters. */
void Topology::AssignBondParams(ParmHolder<BondParmType> const& newBondParams) {
  bondparm_.clear();
  ParmHolder<int> currentIndices;
  //AssignParm<BondParmType, BondArray, BondParmArray>( atoms, newBondParams, currentIndices, bonds_, bondparm_ );
  AssignBondParm( newBondParams, currentIndices, bonds_,  bondparm_, "bond" );
  AssignBondParm( newBondParams, currentIndices, bondsh_, bondparm_, "bond" );
}

/** Replace any current Urey-Bradley parameters with given UB parameters. */
void Topology::AssignUBParams(ParmHolder<BondParmType> const& newBondParams) {
  chamber_.SetUBparm().clear();
  ParmHolder<int> currentIndices;
  AssignBondParm( newBondParams, currentIndices, chamber_.SetUB(), chamber_.SetUBparm(), "UB term" );
}

/** Set parameters for angles in given angle array. */
void Topology::AssignAngleParm(ParmHolder<AngleParmType> const& newAngleParams,
                              ParmHolder<int>& currentIndices,
                              AngleArray& angles)
{
  for (AngleArray::iterator ang = angles.begin(); ang != angles.end(); ++ang) {
    TypeNameHolder types(3);
    types.AddName( atoms_[ang->A1()].Type() );
    types.AddName( atoms_[ang->A2()].Type() );
    types.AddName( atoms_[ang->A3()].Type() );
    bool found;
    // See if parameter already present.
    int idx = currentIndices.FindParam( types, found );
    if (!found) {
      // Search in new
      AngleParmType ap = newAngleParams.FindParam( types, found );
      if (found) {
        // Add parameter
        idx = (int)angleparm_.size();
        angleparm_.push_back( ap );
        currentIndices.AddParm( types, idx, false );
      } else
        idx = -1;
    }
    if (idx == -1)
      mprintf("Warning: Angle parameter not found for angle %s-%s-%s (%s-%s-%s)\n",
              TruncResAtomNameNum(ang->A1()).c_str(),
              TruncResAtomNameNum(ang->A2()).c_str(),
              TruncResAtomNameNum(ang->A3()).c_str(),
              *types[0], *types[1], *types[3]);
    ang->SetIdx( idx );
  }
}

/** Replace any current angle parameters with given angle parameters. */
void Topology::AssignAngleParams(ParmHolder<AngleParmType> const& newAngleParams) {
  angleparm_.clear();
  ParmHolder<int> currentIndices;
  AssignAngleParm( newAngleParams, currentIndices, angles_ );
  AssignAngleParm( newAngleParams, currentIndices, anglesh_ );
}

/** Set parameters for dihedrals in given dihedral array. */
void Topology::AssignImproperParm(ParmHolder<DihedralParmType> const& newImproperParams,
                                  ParmHolder<int>& currentIndices,
                                  DihedralArray& impropers)
{
  for (DihedralArray::iterator imp = impropers.begin(); imp != impropers.end(); ++imp) {
    TypeNameHolder types(4);
    types.AddName( atoms_[imp->A1()].Type() );
    types.AddName( atoms_[imp->A2()].Type() );
    types.AddName( atoms_[imp->A3()].Type() );
    types.AddName( atoms_[imp->A4()].Type() );
    bool found;
    // See if parameter already present.
    int idx = currentIndices.FindParam( types, found );
    if (!found) {
      // Search in new
      DihedralParmType ip = newImproperParams.FindParam( types, found );
      if (found) {
        // Add parameter
        idx = (int)chamber_.ImproperParm().size();
        chamber_.SetImproperParm().push_back( ip );
        currentIndices.AddParm( types, idx, false );
      } else
        idx = -1;
    }
    if (idx == -1)
      mprintf("Warning: Parameter not found for improper %s-%s-%s-%s (%s-%s-%s-%s)\n",
              TruncResAtomNameNum(imp->A1()).c_str(),
              TruncResAtomNameNum(imp->A2()).c_str(),
              TruncResAtomNameNum(imp->A3()).c_str(),
              TruncResAtomNameNum(imp->A4()).c_str(),
              *types[0], *types[1], *types[3], *types[4]);
    imp->SetIdx( idx );
  }
}

/** Replace any current improper parameters with given improper parameters. */
void Topology::AssignImproperParams(ParmHolder<DihedralParmType> const& newImproperParams) {
  chamber_.SetImproperParm().clear();
  ParmHolder<int> currentIndices;
  AssignImproperParm( newImproperParams, currentIndices, chamber_.SetImpropers() );
}

/** Set parameters for dihedrals in given dihedral array. */
void Topology::AssignDihedralParm(DihedralParmHolder const& newDihedralParams,
                                  DihedralArray& dihedralsIn)
{
  // Dihedrals can be a bit of a pain since there can be multiple
  // multiplicities for a single dihedral type. In case multiplicities
  // change, start with a fresh dihedral array containing only unique
  // dihedrals.
  DihedralArray tmpdih = dihedralsIn;
  std::sort( tmpdih.begin(), tmpdih.end() );
  DihedralArray dihedrals;
  for (DihedralArray::iterator dih = tmpdih.begin(); dih != tmpdih.end(); ++dih) {
    if (dihedrals.empty())
      dihedrals.push_back( *dih );
    else {
      DihedralType const& last = dihedrals.back();
      if ( last.A1() != dih->A1() ||
           last.A2() != dih->A2() ||
           last.A3() != dih->A3() ||
           last.A4() != dih->A4() )
        dihedrals.push_back( *dih );
    }
  }
  if (debug_ > 0)
    mprintf("DEBUG: %zu incoming dihedrals, %zu unique dihedrals.\n",
            dihedralsIn.size(), dihedrals.size());

  ParmHolder< std::vector<int> > currentIndices;
  dihedralsIn.clear();
  for (DihedralArray::iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih) {
    TypeNameHolder types(4);
    types.AddName( atoms_[dih->A1()].Type() );
    types.AddName( atoms_[dih->A2()].Type() );
    types.AddName( atoms_[dih->A3()].Type() );
    types.AddName( atoms_[dih->A4()].Type() );
    bool found;
    // See if parameter already present.
    std::vector<int> idxs = currentIndices.FindParam( types, found );
    if (!found) {
      // Not yet present; search in newDihedralParams for parameter(s).
      DihedralParmArray dpa = newDihedralParams.FindParam( types, found );
      if (found) {
        for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
        {
          // Add parameter and save the index
          int idx = (int)dihedralparm_.size();
          idxs.push_back( idx );
          dihedralparm_.push_back( *it );
        }
        currentIndices.AddParm( types, idxs, false );
      }
    }
    if (idxs.empty()) {
      mprintf("Warning: Dihedral parameters not found for dihedral %s-%s-%s-%s (%s-%s-%s-%s)\n",
              TruncResAtomNameNum(dih->A1()).c_str(),
              TruncResAtomNameNum(dih->A2()).c_str(),
              TruncResAtomNameNum(dih->A3()).c_str(),
              TruncResAtomNameNum(dih->A4()).c_str(),
              *types[0], *types[1], *types[3], *types[4]);
      DihedralType mydih = *dih;
      mydih.SetIdx( -1 );
      dihedralsIn.push_back( mydih );
    } else {
      bool multi = idxs.size() > 1;
      // Actually add the dihedrals
      for (std::vector<int>::const_iterator dpidx = idxs.begin(); dpidx != idxs.end(); ++dpidx)
      {
        DihedralType mydih = *dih;
        mydih.SetIdx( *dpidx );
        if (multi) {
          // If there are multiple parameters for the same dihedral, all but
          // one of the 1-4 calcs need to be skipped.
          if (dpidx+1 != idxs.end())
            mydih.SetSkip14(true);
          else
            mydih.SetSkip14(false);
        }
        dihedralsIn.push_back( mydih );
      }
    }
  }
}

/** Replace any current dihedral parameters with given dihedral parameters. */
void Topology::AssignDihedralParams(DihedralParmHolder const& newDihedralParams) {
  dihedralparm_.clear();
  AssignDihedralParm( newDihedralParams, dihedrals_ );
  AssignDihedralParm( newDihedralParams, dihedralsh_ );
}

/** Replace current nonbond parameters with given nonbond parameters. */
//TODO Accept array of atom type names?
void Topology::AssignNonbondParams(ParmHolder<AtomType> const& newTypes, ParmHolder<NonbondType> const& newNB) {
  // Generate array of existing types.
  ParmHolder<AtomType> atomTypes;
  for (std::vector<Atom>::const_iterator atm = atoms_.begin(); atm != atoms_.end(); ++atm)
  {
    if (atm->Type().len() > 0) {
      TypeNameHolder types(1);
      types.AddName(atm->Type());
      // Find in newTypes.
      bool found;
      AtomType newAT = newTypes.FindParam( types, found );
      if (!found) {
        mprintf("Warning: No atom type information for type '%s'\n", *types[0]);
        newAT = AtomType(atm->Mass());
      }
      atomTypes.AddParm( types, newAT, true );
    }
  }
  if (atomTypes.size() < 1) {
    mprintf("Warning: No atom type information in %s - cannot assign nonbond parameters.\n",
            c_str());
    return;
  }
  // Regenerate nonbond params for existing types
  nonbond_.Clear();
  nonbond_.SetupLJforNtypes( atomTypes.size() );
  int nidx1 = 0;
  NameIdxMapType nameIdxMap;
  for (ParmHolder<AtomType>::const_iterator t1 = atomTypes.begin(); t1 != atomTypes.end(); ++t1, nidx1++)
  {
    NameType const& name1 = t1->first[0];
    //mprintf("DEBUG: Type1= %s (%i)\n", *name1, nidx1);
    AtomType const& type1 = t1->second;
    nameIdxMap.insert( std::pair<NameType, int>(name1, nidx1) );
    int nidx2 = nidx1;
    for (ParmHolder<AtomType>::const_iterator t2 = t1; t2 != atomTypes.end(); ++t2, nidx2++)
    {
      NameType const& name2 = t2->first[0];
      //mprintf("DEBUG:\t\tType2= %s (%i)\n", *name2, nidx2);
      AtomType const& type2 = t2->second;
      TypeNameHolder types(2);
      types.AddName( name1 );
      types.AddName( name2 );
      // See if this parameter exists in the given nonbond array.
      NonbondType LJAB;
      ParmHolder<NonbondType>::const_iterator it = newNB.GetParam( types );
      if (it == newNB.end()) {
        mprintf("NB parameter for %s %s not found. Generating.\n", *name1, *name2);
        LJAB = type1.LJ().Combine_LB( type2.LJ() );
      } else {
        mprintf("Using existing NB parameter for %s %s\n", *name1, *name2);
        LJAB = it->second;
      }
      nonbond_.AddLJterm(nidx1, nidx2, LJAB);
    }
  }
  // Reset the atom type indices.
  for (std::vector<Atom>::iterator atm = atoms_.begin(); atm != atoms_.end(); ++atm)
  {
    int tidx = -1;
    NameIdxMapType::const_iterator it = nameIdxMap.find( atm->Type() );
    if (it == nameIdxMap.end()) {
      mprintf("Warning: Atom type not found for %s (type %s)\n",
              TruncResAtomNameNum( atm - atoms_.begin() ).c_str(),
              *(atm->Type()));
    } else {
      if (it->second < 0 || it->second >= (int)atomTypes.size()) {
        mprinterr("Internal Error: Type index for %s (type %s) out of range: %i\n",
                  TruncResAtomNameNum( atm - atoms_.begin() ).c_str(),
                  *(atm->Type()), it->second);
      } else
        tidx = it->second;
    }
    atm->SetTypeIndex( tidx );
  }
}

// -----------------------------------------------------------------------------

/** Update/add to parameters in this topology with those from given set. */
int Topology::UpdateParams(ParameterSet const& set1) {
  ParameterSet set0 = GetParameters();
  // Check
  if (set0.AT().size() < 1)
    mprintf("Warning: No atom type information in '%s'\n", c_str());
  if (debug_ > 0) {
    mprintf("DEBUG: Saving original parameters in originalp.dat, new parameters in newp.dat\n");
    set0.Debug("originalp.dat");
  }

  unsigned int updateCount;
  // Bond parameters
  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.BP(), set1.BP(), "bond");
  if (updateCount > 0) {
    mprintf("\tRegenerating bond parameters.\n");
    AssignBondParams( set0.BP() );
  }
  // Angle parameters
  updateCount = UpdateParameters< ParmHolder<AngleParmType> >(set0.AP(), set1.AP(), "angle");
  if (updateCount > 0) {
    mprintf("\tRegenerating angle parameters.\n");
    AssignAngleParams( set0.AP() );
  }
  // Dihedral parameters
  updateCount = UpdateParameters< DihedralParmHolder >(set0.DP(), set1.DP(), "dihedral");
  if (updateCount > 0) {
    mprintf("\tRegenerating dihedral parameters.\n");
    AssignDihedralParams( set0.DP() );
  }
  // Urey-Bradley
  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.UB(), set1.UB(), "Urey-Bradley");
  if (updateCount > 0) {
    mprintf("\tRegenerating UB parameters.\n");
    AssignUBParams( set0.UB() );
  }
  // Improper parameters
  updateCount = UpdateParameters< ParmHolder<DihedralParmType> >(set0.IP(), set1.IP(), "improper");
  if (updateCount > 0) {
    mprintf("\tRegenerating improper parameters.\n");
    AssignImproperParams( set0.IP() );
  }
  // Atom types
  updateCount = UpdateParameters< ParmHolder<AtomType> >(set0.AT(), set1.AT(), "atom type");
  updateCount += UpdateParameters< ParmHolder<NonbondType> >(set0.NB(), set1.NB(), "LJ A-B");
  if (updateCount > 0) {
    mprintf("\tRegenerating nonbond parameters.\n");
    AssignNonbondParams( set0.AT(), set0.NB() );
  }

  if (debug_ > 0) set0.Debug("newp.dat");
  return 0;
}
