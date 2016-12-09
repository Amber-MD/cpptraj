#include <cmath> // pow
#include <algorithm> // find
#include <stack> // For large system molecule search
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
#include "Constants.h" // RADDEG, SMALL

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
  * "<resname><resnum>@<atomname>", e.g. "ARG_11@CA".
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

// Topology::PrintDihedrals()
void Topology::PrintDihedrals(DihedralArray const& darray, CharMask const& maskIn, 
                              int& nd, bool Select_OR) const
{
  int rwidth = DigitWidth(residues_.size()) + 7;
  for (DihedralArray::const_iterator datom = darray.begin();
                                     datom != darray.end(); ++datom)
  {
    int atom1 = (*datom).A1();
    int atom2 = (*datom).A2();
    int atom3 = (*datom).A3();
    int atom4 = (*datom).A4();
    bool selected;
    if (Select_OR)
      selected = (maskIn.AtomInCharMask( atom1 ) || maskIn.AtomInCharMask( atom2 ) ||
                  maskIn.AtomInCharMask( atom3 ) || maskIn.AtomInCharMask( atom4 )   );
    else // AND
      selected = (maskIn.AtomInCharMask( atom1 ) && maskIn.AtomInCharMask( atom2 ) &&
                  maskIn.AtomInCharMask( atom3 ) && maskIn.AtomInCharMask( atom4 )   );
    if (selected) {
      // Determine dihedral type: 'E'nd, 'I'mproper, or 'B'oth
      char type = ' ';
      if ((*datom).Type() == DihedralType::END) type = 'E';
      else if ((*datom).Type() == DihedralType::IMPROPER) type = 'I';
      else if ((*datom).Type() == DihedralType::BOTH) type = 'B';
      mprintf("%c %8i:", type, nd);
      int didx = (*datom).Idx();
      if ( didx > -1 )
        mprintf(" %6.3f %4.2f %4.1f", dihedralparm_[didx].Pk(), dihedralparm_[didx].Phase(),
                 dihedralparm_[didx].Pn());
      mprintf(" %-*s %-*s %-*s %-*s (%i,%i,%i,%i)",
              rwidth, AtomMaskName(atom1).c_str(), rwidth, AtomMaskName(atom2).c_str(), 
              rwidth, AtomMaskName(atom3).c_str(), rwidth, AtomMaskName(atom4).c_str(),
              atom1+1, atom2+1, atom3+1, atom4+1);
      // Atom types
      const char* atype1 = *atoms_[atom1].Type();
      const char* atype2 = *atoms_[atom2].Type();
      const char* atype3 = *atoms_[atom3].Type();
      const char* atype4 = *atoms_[atom4].Type();
      mprintf(" %c%c-%c%c-%c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1],
              atype3[0],atype3[1],atype4[0],atype4[1]);
    }
    nd++;
  }
  mprintf("\n");
}

// Topology::PrintDihedralInfo()
void Topology::PrintDihedralInfo(std::string const& maskString, bool select_OR) const {
  CharMask mask( maskString );
  if (SetupCharMask( mask )) return;
  mprintf("#");
  mask.MaskInfo();
  if (mask.None()) return;
  mprintf("#Dihedral    pk     phase pn                atoms\n");
  int nd = 1;
  if (!dihedralsh_.empty())
    PrintDihedrals( dihedralsh_, mask, nd, select_OR );
  if (!dihedrals_.empty())
    PrintDihedrals( dihedrals_, mask, nd, select_OR );
}


// Topology::PrintMoleculeInfo()
void Topology::PrintMoleculeInfo(std::string const& maskString) const {
  if (molecules_.empty())
    mprintf("\t'%s' No molecule info.\n",c_str());
  else {
    CharMask mask( maskString );
    if (SetupCharMask( mask )) return;
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      int mwidth = DigitWidth(molecules_.size());
      if (mwidth < 5) mwidth = 5;
      int awidth = DigitWidth(atoms_.size());
      if (awidth < 5) awidth = 5;
      int rwidth = DigitWidth(residues_.size());
      if (rwidth < 5) rwidth = 5;
      mprintf("%-*s %*s %*s %4s\n", mwidth, "#Mol", awidth, "Natom", 
              rwidth, "#Res", "Name");
      unsigned int mnum = 1;
      for (std::vector<Molecule>::const_iterator mol = molecules_.begin(); 
                                                 mol != molecules_.end(); mol++)
      {
        if ( mask.AtomsInCharMask( mol->BeginAtom(), mol->EndAtom() ) ) {
          int firstres = atoms_[ mol->BeginAtom() ].ResNum();
          mprintf("%*u %*i %*i %4s %c", mwidth, mnum, awidth, mol->NumAtoms(),
                  rwidth, firstres+1, residues_[firstres].c_str(), residues_[firstres].ChainID());
          if ( mol->IsSolvent() ) mprintf(" SOLVENT");
          mprintf("\n");
        }
        ++mnum;
      }
    }
  }
}

// Topology::PrintResidueInfo()
/** Since this function may be called from command line with worldsilent
  * set to true, use loudPrintf and mprinterr.
  */
void Topology::PrintResidueInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  if (SetupIntegerMask( mask )) return;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    int awidth = DigitWidth(atoms_.size());
    if (awidth < 5) awidth = 5;
    int rwidth = DigitWidth(residues_.size());
    if (rwidth < 5) rwidth = 5;
    int mwidth = DigitWidth(molecules_.size());
    if (mwidth < 5) mwidth = 5;
    loudPrintf("%-*s %4s %*s %*s %*s %*s %*s\n", rwidth, "#Res", "Name",
               awidth, "First", awidth, "Last", 
               awidth, "Natom", rwidth, "#Orig", mwidth, "#Mol");
    int rn = -1;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      if (atoms_[*atom].ResNum() > rn) {
        rn = atoms_[*atom].ResNum();
        Residue const& res = residues_[rn];
        loudPrintf("%*i %4s %*i %*i %*i %*i %*i %c\n", rwidth, rn+1, res.c_str(),
                   awidth, res.FirstAtom()+1, awidth, res.LastAtom(),
                   awidth, res.NumAtoms(), rwidth, res.OriginalResNum(),
                   mwidth, atoms_[*atom].MolNum()+1, res.ChainID());
      }
    }
  }
}

/** Print residue info using single char names. */
void Topology::PrintShortResInfo(std::string const& maskString, int maxChar) const {
  AtomMask mask( maskString );
  if (SetupIntegerMask( mask )) return;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    // Determine last selected residue.
    int max_res = atoms_[mask.back()].ResNum();
    int total = 0;
    int rn = -1, startRes = -1;
    std::string resLine;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      int current_res = atoms_[*atom].ResNum();
      if (current_res > rn) {
        int n_res_skipped = 1;
        if (startRes == -1)
          startRes = current_res;
        else
          n_res_skipped = current_res - rn;
        // If we skipped any residues print last consective segment and start a new one.
        if (n_res_skipped > 1) {
          mprintf("%-8i %s\n", startRes+1, resLine.c_str());
          startRes = current_res;
          resLine = residues_[current_res].SingleCharName();
          total = 1;
        } else {
          // Convert residue name.
          resLine += residues_[current_res].SingleCharName();
          total++;
        }
        // Print if max line length reached or final res.
        if ((total%maxChar)==0 || current_res == max_res)
        {
          mprintf("%-8i %s\n", startRes+1, resLine.c_str());
          if (current_res == max_res) break;
          startRes = -1;
          resLine.clear();
        } else if ((total % 10) == 0)
          resLine += ' ';
        rn = current_res;
      }
    }
  }
}

// Topology::PrintChargeMassInfo()
int Topology::PrintChargeMassInfo(std::string const& maskString, int type) const {
  AtomMask mask( maskString );
  if (SetupIntegerMask( mask )) return 1;
  if (type == 0 || type == 2) {
    mprintf("\tSum of charges in mask");
    mask.BriefMaskInfo();
    double sumq = 0.0;
    for (AtomMask::const_iterator aidx = mask.begin(); aidx != mask.end(); ++aidx)
      sumq += atoms_[*aidx].Charge();
    mprintf(" is %g\n", sumq);
  }
  if (type == 1 || type == 2) {
    mprintf("\tSum of masses in mask");
    mask.BriefMaskInfo();
    double summ = 0.0;
    for (AtomMask::const_iterator aidx = mask.begin(); aidx != mask.end(); ++aidx)
      summ += atoms_[*aidx].Mass();
    mprintf(" is %g\n", summ);
  }
  return 0; 
}

// -----------------------------------------------------------------------------
// Topology::AddTopAtom()
int Topology::AddTopAtom(Atom const& atomIn, Residue const& resIn)
{
  // If no residues or res num has changed, this is a new residue.
  if ( residues_.empty() || 
       residues_.back().OriginalResNum() != resIn.OriginalResNum() ||
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
          break;
        }
      }
    }
    if (mols_share_residues) {
      mprintf("Warning: 2 or more molecules share residue numbers.\n");
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
/** Create a bond between atom1 and atom2, update the atoms array.
  * For bonds to H always insert the H second.
  */
void Topology::AddBond(int atom1, int atom2, int pidxIn) {
  // Check if atoms are out of range.
  if (atom1 < 0 || atom1 >= (int)atoms_.size()) {
    mprintf("Warning: Atom # %i is out of range, cannot create bond.\n", atom1+1);
    return;
  }
  if (atom2 < 0 || atom2 >= (int)atoms_.size()) {
    mprintf("Warning: Atom # %i is out of range, cannot create bond.\n", atom2+1);
    return;
  }
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

void Topology::AddAngle(int atom1, int atom2, int atom3) {
  // FIXME: Check duplicate
  if (atoms_[atom1].Element() == Atom::HYDROGEN ||
      atoms_[atom2].Element() == Atom::HYDROGEN ||
      atoms_[atom3].Element() == Atom::HYDROGEN)
    anglesh_.push_back( AngleType(atom1, atom2, atom3, -1) );
  else
    angles_.push_back( AngleType(atom1, atom2, atom3, -1) );
}

void Topology::AddAngle(AngleType const& angIn, bool isH) {
  if (isH)
    anglesh_.push_back( angIn );
  else
    angles_.push_back( angIn );
}

void Topology::AddDihedral(int atom1, int atom2, int atom3, int atom4) {
  // FIXME: Check duplicate
  if (atoms_[atom1].Element() == Atom::HYDROGEN ||
      atoms_[atom2].Element() == Atom::HYDROGEN ||
      atoms_[atom3].Element() == Atom::HYDROGEN ||
      atoms_[atom4].Element() == Atom::HYDROGEN)
    dihedralsh_.push_back( DihedralType(atom1, atom2, atom3, atom4, -1) );
  else
    dihedrals_.push_back( DihedralType(atom1, atom2, atom3, atom4, -1) );
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
  return mask.SetupMask(atoms_, residues_, refCoords_.xAddress());
}

// Topology::SetupCharMask()
int Topology::SetupCharMask(CharMask &mask) const {
  return mask.SetupMask(atoms_, residues_, refCoords_.xAddress());
}

// Topology::SetupIntegerMask()
int Topology::SetupIntegerMask(AtomMask &mask, Frame const& frame) const {
  if (frame.empty()) return mask.SetupMask(atoms_, residues_, 0);
  return mask.SetupMask(atoms_, residues_, frame.xAddress());
}

// Topology::SetupCharMask()
int Topology::SetupCharMask(CharMask &mask, Frame const& frame) const {
  if (frame.empty()) return mask.SetupMask(atoms_, residues_, 0);
  return mask.SetupMask(atoms_, residues_, frame.xAddress());
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
  // Set new solvent information based on new molecules
  if (newParm->SetSolventInfo()) {
    mprintf("Warning: Could not set up solvent information for stripped topology %s\n",
            newParm->c_str());
  } 
  // Set up new angle info
  newParm->angles_ = StripAngleArray( angles_, atomMap );
  newParm->anglesh_ = StripAngleArray( anglesh_, atomMap );
  parmMap.assign( angleparm_.size(), -1 );
  StripAngleParmArray( newParm->angles_,  parmMap, newParm->angleparm_ );
  StripAngleParmArray( newParm->anglesh_, parmMap, newParm->angleparm_ );
  // Set up new dihedral info
  newParm->dihedrals_ = StripDihedralArray( dihedrals_, atomMap );
  newParm->dihedralsh_ = StripDihedralArray( dihedralsh_, atomMap );
  parmMap.assign( dihedralparm_.size(), -1 );
  StripDihedralParmArray( newParm->dihedrals_,  parmMap, newParm->dihedralparm_ );
  StripDihedralParmArray( newParm->dihedralsh_, parmMap, newParm->dihedralparm_ );
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
    for (int a1idx = 0; a1idx != (int)oldTypeArray.size(); a1idx++)
    {
      int atm1 = oldTypeArray[a1idx];
      for (int a2idx = a1idx; a2idx != (int)oldTypeArray.size(); a2idx++)
      {
        int atm2 = oldTypeArray[a2idx];
        int oldnbidx = nonbond_.GetLJindex( atm1, atm2 );
        // NOTE: Certain routines in sander (like the 1-4 calcs) do NOT use
        //       the nonbond index array; instead they expect the nonbond
        //       arrays to be indexed like '(ibig*(ibig-1)/2+isml)', where
        //       ibig is the larger atom type index.
        int ibig = std::max(a1idx, a2idx) + 1;
        int isml = std::min(a1idx, a2idx) + 1;
        int testidx = (ibig*(ibig-1)/2+isml)-1;
        if (oldnbidx > -1) {
          // This is a traditional LJ 6-12 term. Because of the way the LJ 1-4
          // code is laid out in sander/pmemd the LJ matrix has to be laid out
          // indepdendent of the nonbond index array.
          newParm->nonbond_.AddLJterm( testidx, a1idx, a2idx, nonbond_.NBarray(oldnbidx) );
        } else {
          // This is an old LJ 10-12 hbond term. Add one to the LJ 6-12 matrix
          // and one to the hbond since that seems to be the convention.
          newParm->nonbond_.AddLJterm( testidx, a1idx, a2idx, NonbondType() );
          newParm->nonbond_.AddHBterm( a1idx, a2idx, nonbond_.HBarray((-oldnbidx)-1) );
        }
        //int newnbidx = newParm->nonbond_.GetLJindex( a1idx, a2idx );
        //mprintf("DEBUG: oldtypei=%i oldtypej=%i Old NB index=%i, newtypi=%i newtypej=%i new NB idx=%i testidx=%i\n", 
        //        atm1, atm2, oldnbidx, a1idx, a2idx, newnbidx, testidx);
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
  // CHAMBER info - Parameters remain intact
  if (chamber_.HasChamber()) {
    newParm->chamber_.SetVersion( chamber_.FF_Version(), chamber_.FF_Type() );
    newParm->chamber_.SetUB( StripBondArray(chamber_.UB(),atomMap), chamber_.UBparm() );
    newParm->chamber_.SetImproper( StripDihedralArray(chamber_.Impropers(),atomMap),
                                   chamber_.ImproperParm() );
    newParm->chamber_.SetLJ14( chamber_.LJ14() );
    if (chamber_.HasCmap()) {
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
  for (BondArray::iterator bnd = newBondArray.begin();
                           bnd != newBondArray.end(); ++bnd)
  {
    int oldidx = bnd->Idx();
    int newidx = parmMap[bnd->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newBondParm.size();
      parmMap[oldidx] = newidx;
      newBondParm.push_back( bondparm_[oldidx] );
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
  for (DihedralArray::iterator dih = newDihedralArray.begin();
                               dih != newDihedralArray.end(); ++dih)
  {
    int oldidx = dih->Idx();
    int newidx = parmMap[dih->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newDihedralParm.size();
      parmMap[oldidx] = newidx;
      newDihedralParm.push_back( dihedralparm_[oldidx] );
    }
    //mprintf("DEBUG: Old dihedral parm index=%i, new dihedral parm index=%i\n", oldidx, newidx);
    dih->SetIdx( newidx );
  }
}

// Topology::AddBondArray()
void Topology::AddBondArray(BondArray const& barray, int atomOffset) {
  for (BondArray::const_iterator bond = barray.begin(); bond != barray.end(); ++bond)
    AddBond( bond->A1() + atomOffset, bond->A2() + atomOffset );
}

// Topology::AppendTop()
int Topology::AppendTop(Topology const& CurrentTop) {
  int atomOffset = (int)atoms_.size();
  int resOffset = (int)residues_.size();
  // ATOMS
  for (atom_iterator atom = CurrentTop.begin(); atom != CurrentTop.end(); ++atom)
  {
    Atom CurrentAtom = *atom;
    Residue const& res = CurrentTop.Res( CurrentAtom.ResNum() );
    // Bonds need to be cleared and re-added.
    CurrentAtom.ClearBonds();
    AddTopAtom( CurrentAtom, Residue(res.Name(), CurrentAtom.ResNum() + resOffset,
                                     res.Icode(), res.ChainID()) );
  }
  // BONDS
  AddBondArray(CurrentTop.Bonds(),  atomOffset);
  AddBondArray(CurrentTop.BondsH(), atomOffset);
  // Re-set up this topology
  // TODO: Could get expensive for multiple appends.
  return CommonSetup();
}
