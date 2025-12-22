#include <algorithm> // find, copy, fill
#include <stack> // For large system molecule search
#include "Topology.h"
#include "AtomMask.h"
#include "AtomType.h" // AppendTop()
#include "CharMask.h"
#include "Constants.h" // SMALL
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "Parm/AssignParams.h" // AppendTop()
#include "Parm/GetParams.h" // AppendTop()
#include "Parm/Merge.h" // AppendTop()
#include "Parm/ParmHolder.h" // AppendTop()
#include "UpdateParameters.h" // AppendTop(). This include must be last

const NonbondType Topology::LJ_EMPTY = NonbondType();

// CONSTRUCTOR
Topology::Topology() :
  debug_(0),
  ipol_(0),
  NsolventMolecules_(0),
  pindex_(0),
  n_extra_pts_(0)
{ }

/** Copy the metadata from another Topology. */
void Topology::CopyTopMetadata(Topology const& rhs) {
  parmName_ = rhs.parmName_;
  fileName_ = rhs.fileName_;
  ff_desc_ = rhs.ff_desc_;
  parmBox_ = rhs.parmBox_;
}

/** Set the parm name only. */
void Topology::SetParmTitle(std::string const& title) {
  parmName_ = title;
}

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

/** \return count of heavy atoms (ignore hydrogen/extra points). */
unsigned int Topology::HeavyAtomCount() const {
  unsigned int hac = 0;
  for (atom_iterator at = begin(); at != end(); ++at) {
    if (at->Element() != Atom::HYDROGEN &&
        at->Element() != Atom::EXTRAPT)
      hac++;
  }
  return hac;
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
    res->SetChainID("");
  }
  missingRes_.clear();
  missingHet_.clear();
}

/** Set list of missing residues and residues missing heteroatoms */
void Topology::SetMissingResInfo(std::vector<Residue> const& missingResIn,
                                 std::vector<Residue> const& missingHetIn)
{
  missingRes_ = missingResIn;
  missingHet_ = missingHetIn;
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

/** Merge consecutive residues into a single residue. */
int Topology::MergeResidues(int startRes, int stopRes) {
  // Check that start and stop make sense
  if (startRes < 0 || stopRes < 0) {
    mprinterr("Internal Error: MergeResidues: Either startRes (%i) or stopRes (%i) < 0.\n",
              startRes, stopRes);
    return 1;
  }
  if (stopRes < startRes) {
    mprinterr("Error: Start residue %i > stop residue %i. Cannot merge residues.\n",
              startRes+1, stopRes+1);
    return 1;
  } else if (startRes == stopRes) {
    mprintf("Warning: Start residue %i is the same as stop residue %i. Nothing to merge.\n",
            startRes+1, stopRes+1);
    return 0;
  }
  // TODO check that residues are bonded?
  // Check for duplicate atom names
  int startAtom = Res(startRes).FirstAtom();
  int stopAtom = Res(stopRes).LastAtom();
  for (int at0 = startAtom; at0 != stopAtom; at0++) {
    for (int at1 = at0+1; at1 != stopAtom; at1++) {
      if ( (*this)[at0].Name() == (*this)[at1].Name() ) {
        mprintf("Warning: In merge of residues %i-%i, duplicate atom name %s (%s)\n",
                startRes+1, stopRes+1,
                AtomMaskName(at0).c_str(), AtomMaskName(at1).c_str());
      }
    }
  }
  std::vector<Residue> newResidues;
  unsigned int newNres = residues_.size() - (unsigned int)(stopRes - startRes);
  mprintf("\tOriginally %zu residues, merging %i-%i, now %u residues.\n",
          residues_.size(), startRes+1, stopRes+1, newNres);
  newResidues.reserve( newNres );
  // Residues up to but not including startRes
  for (int ires = 0; ires != startRes; ires++)
    newResidues.push_back( Res(ires) );
  // The merged residue. Use chain ID etc of the first residue.
  Residue mergedRes = Res(startRes);
  mergedRes.SetLastAtom( Res(stopRes).LastAtom() );
  newResidues.push_back( mergedRes );
  // Update the atoms of the merged residues
  for (int at = Res(startRes+1).FirstAtom(); at != Res(stopRes).LastAtom(); at++)
    atoms_[at].SetResNum( startRes );
  // Residues from after stopRes to end
  for (int ires = stopRes + 1; ires < Nres(); ires++)
    newResidues.push_back( Res(ires) );
  // Overwrite old residue info
  residues_ = newResidues;

  return 0;
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

/** Given an atom number, return a string in LEaP-style format. */
std::string Topology::LeapName(int at) const {
  int rnum = atoms_[at].ResNum();
  int idx = at - residues_[rnum].FirstAtom() + 1;
  std::string out( ".R<" + residues_[rnum].Name().Truncated() + " " + integerToString(rnum+1) +
                  ">.A<" + atoms_[at].Name().Truncated()      + " " + integerToString(idx) + ">");
  return out;
}

// Topology::TruncResNameNum()
/** Given a residue index (starting from 0), return a string containing 
  * residue name and number (starting from 1) with format: 
  * "<resname>:<resnum>", e.g. "ARG:11".
  * Truncate residue name so there are no blanks.
  */
std::string Topology::TruncResNameNum(int res) const {
  if (res < 0 || res >= (int)residues_.size()) return std::string("");
  // Residue name with no trailing spaces.
  return residues_[res].Name().Truncated() + ":" + integerToString( res+1 );
}

/** Given a residue index, return a string containing residue name,
  * original residue number, and (optionally) chain ID with format:
  * "<resname>_<onum>[_<id>]".
  * Truncate residue name so there are no blanks.
  */
std::string Topology::TruncResNameOnumId(int res) const {
  if (res < 0 || res >= (int)residues_.size()) return std::string("");
  std::string name = residues_[res].Name().Truncated() + "_" +
                     integerToString(residues_[res].OriginalResNum());
  if (residues_[res].HasChainID())
    name.append( "_" + residues_[res].ChainID() );
  return name;
}

/** Given an atom index, return a string containing the residue name,
  * original residue number, an optional chain ID, and the atom
  * name with format:
  * "<resname>_<onum>[_<id>]@<atom name>".
  * Truncate residue and atom names so there are no blanks.
  */
std::string Topology::TruncAtomResNameOnumId(int at) const {
  if (at < 0 || at >= (int)atoms_.size()) return std::string("");
  Atom const& thisAtom = atoms_[at];
  std::string name = TruncResNameOnumId( thisAtom.ResNum() );
  name.append( "@" + thisAtom.Name().Truncated() );
  return name;
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
  if (HasChamber()) {
    mprintf("\t\tCHAMBER: %zu Urey-Bradley terms, %zu Impropers, %zu LJ 1-4 terms.\n",
            ub_.size(), impropers_.size(), nonbond_.LJ14().size());
  }
  if (HasCmap())
    mprintf("\t\t%zu CMAP grids, %zu CMAP terms.\n", 
            CmapGrid().size(), Cmap().size());
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

/** This version takes a molecule number as well. */
int Topology::addTopAtom(Atom const& atomIn, Residue const& resIn,
                         unsigned int molnum, bool isSolvent)
{
  unsigned int atnum = atoms_.size();
  AddTopAtom(atomIn, resIn);
  atoms_.back().SetMol( molnum );
  if (molnum > molecules_.size())
    mprintf("Warning: AddTopAtom(): Molecule number %i is not consecutive.\n");
  if (molnum >= molecules_.size()) {
    molecules_.resize( molnum+1 );
    if (isSolvent)
      NsolventMolecules_++;
  }
  molecules_[molnum].AddAtnum( atnum );
  return 0;
}

/** Add specified solvent residues from given topology to this topology. */
int Topology::AddSolventResidues(Topology const& solventTop, std::vector<int> const& solventResNums,
                                 Frame& frameOut, Frame const& solventFrame)
{
  std::vector<int> bondedAtoms;
  for (std::vector<int>::const_iterator ires = solventResNums.begin();
                                        ires != solventResNums.end(); ++ires)
  {
    int atomOffset = Natom();
    Residue solventRes = solventTop.Res(*ires);
    solventRes.SetOriginalNum( residues_.size() + 1 );
    bondedAtoms.clear();
    for (int iat = solventRes.FirstAtom(); iat != solventRes.LastAtom(); iat++)
    {
      Atom solventAtom = solventTop[iat];
      // Save solvent bonds
      for (Atom::bond_iterator bat = solventAtom.bondbegin(); bat != solventAtom.bondend(); ++bat) {
        if (*bat > iat) {
          bondedAtoms.push_back(  iat + atomOffset - solventRes.FirstAtom() );
          bondedAtoms.push_back( *bat + atomOffset - solventRes.FirstAtom() );
        }
      }
      solventAtom.ClearBonds(); // FIXME AddTopAtom should clear
      AddTopAtom( solventAtom, solventRes );
      if (solventAtom.Element() == Atom::EXTRAPT)
        n_extra_pts_++;
      // Add PDB info if the topology already has it. FIXME read from incoming top
      if (!bfactor_.empty()) bfactor_.push_back( 0 );
      if (!occupancy_.empty()) occupancy_.push_back( 1 );
      if (!pdbSerialNum_.empty()) pdbSerialNum_.push_back( Natom() );
      if (!atom_altloc_.empty()) atom_altloc_.push_back( ' ' );
      // Add extra arrays
      if (!tree_.empty()) tree_.push_back("BLA");
      if (!ijoin_.empty()) ijoin_.push_back(0);
      if (!irotat_.empty()) irotat_.push_back(0);
      const double* VXYZ = solventFrame.XYZ(iat);
      frameOut.AddXYZ( VXYZ );
    } // END loop over solvent atoms
    // Add bonds
    for (std::vector<int>::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it) {
      int at0 = *it;
      ++it;
      AddBond( at0, *it );
    }
    // Add molecule
    molecules_.push_back( Molecule(atomOffset, Natom()) );
    molecules_.back().SetSolvent();
    NsolventMolecules_++;
  } // END loop over solvent unit residues
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
/** Set up common to all topologies. Assign bond lengths if none are present.
  * \param molsearch If true, determine molecules based on bond info.
  * \param renumberResidues If true, renumber residues if any residue is part of more than 1 molecule
  *        e.g. when alternate locations are present.
  */
int Topology::CommonSetup(bool molsearch, bool renumberResidues)
{
  return CommonSetup(molsearch, renumberResidues, true);
}

/** Set up common to all topologies.
  * \param molsearch If true, determine molecules based on bond info.
  * \param renumberResidues If true, renumber residues if any residue is part of more than 1 molecule
  *        e.g. when alternate locations are present.
  * \param assignBondParm If true, assign default bond lengths if no parameters present.
  */
int Topology::CommonSetup(bool molsearch, bool renumberResidues, bool assignBondParm)
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
  if (assignBondParm && bondparm_.empty())
    generateBondParameters();
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
    residues_.push_back( Residue(res_name, resnum+1, ' ', "") );
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
  ff_desc_.clear();
  ub_.clear();
  ubparm_.clear();
  impropers_.clear();
  improperparm_.clear();
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
/** Check if given bond parm exists in given bond parm array. Add if not.
  * \return Index in bond parm array.
  */
int Topology::addBondParm(BondParmArray& bparray, BondParmType const& BPin)
{
  // See if the BondParm exists.
  int pidx = -1;
  for (BondParmArray::const_iterator bp = bparray.begin();
                                     bp != bparray.end(); ++bp)
  {
    if (BPin == *bp) {
      pidx = (int)(bp - bparray.begin());
      break;
    }
  }
  if (pidx == -1) {
    pidx = (int)bparray.size();
    bparray.push_back( BPin );
  }
  return pidx;
}

/** Create parameters for given bond based on element types. */
void Topology::genBondParam(BondType& bnd, BP_mapType& bpMap)
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

/** Fill in bond parameters based on atomic element types. */
void Topology::generateBondParameters() {
  mprintf("Warning: Determining bond length parameters from element types for '%s'.\n", c_str());
  bondparm_.clear();
  // Hold indices into bondparm for unique element pairs
  BP_mapType bpMap;
  for (BondArray::iterator bnd = bondsh_.begin(); bnd != bondsh_.end(); ++bnd)
    genBondParam( *bnd, bpMap ); 
  for (BondArray::iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
    genBondParam( *bnd, bpMap );
} 

// Topology::AddBond()
void Topology::AddBond(int atom1, int atom2, BondParmType const& BPin) {
  // See if the BondParm exists.
  int pidx = addBondParm( bondparm_, BPin );;
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
  //mprintf("DEBUG: Enter AddBond(%i, %i, %s, %s, %i)\n", atom1+1, atom2+1,
  //        AtomMaskName(atom1).c_str(), AtomMaskName(atom2).c_str(), pidxIn);
  // Check if atoms are out of range.
  if (WarnOutOfRange(atoms_.size(), atom1, "bond")) return;
  if (WarnOutOfRange(atoms_.size(), atom2, "bond")) return;
  bool a1H = (atoms_[atom1].Element() == Atom::HYDROGEN);
  bool a2H = (atoms_[atom2].Element() == Atom::HYDROGEN);
  // Check for duplicate bond
  for (Atom::bond_iterator ba = atoms_[atom1].bondbegin();
                           ba != atoms_[atom1].bondend(); ++ba)
  {
    if ( *ba == atom2 ) {
      if (debug_ > 0)
        mprintf("Warning: Bond between atoms %i and %i already exists.\n", atom1+1, atom2+1);
      // If the bond already exists but does not yet have a parameter, update
      // the parameter index before exiting.
      if (pidxIn > -1 && pidxIn < (int)bondparm_.size()) {
        BondArray* BONDS = 0;
        if (a1H || a2H)
          BONDS = &bondsh_;
        else
          BONDS= &bonds_;
        for (BondArray::iterator it = BONDS->begin(); it != BONDS->end(); ++it) {
          if ( (it->A1() == atom1 && it->A2() == atom2) ||
               (it->A1() == atom2 && it->A2() == atom1) )
          {
            if (debug_ > 0)
              mprintf("DEBUG: Existing bond found. Existing Idx %i Rk=%f Req=%f\n",
                      it->Idx(), bondparm_[it->Idx()].Rk(), bondparm_[it->Idx()].Req());
            if (it->Idx() < 0) {
              if (debug_ > 0)
                mprintf("DEBUG: Adding bond parameter index %i Rk=%f Req=%f for existing bond.\n",
                        pidxIn, bondparm_[pidxIn].Rk(), bondparm_[pidxIn].Req());
              it->SetIdx( pidxIn );
            }
            break;
          }
        } // END loop over target bond array
      }
      return;
    }
  } // END check for duplicate bond.
  // Check if parm index is out of range;
  int pidx;
  if (pidxIn < (int)bondparm_.size())
    pidx = pidxIn;
  else {
    mprintf("Warning: No bond parameters for index %i\n", pidxIn);
    pidx = -1;
  }
  //mprintf("\t\t\tAdding bond %i to %i (a1H=%i a2H=%i) idx=%i\n",atom1+1,atom2+1,(int)a1H,(int)a2H,pidx);
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

/** Clear bond arrays, but do not clear atom connectivity. Used
  * when regenerating bond information from atom connectivity.
  */
void Topology::ClearBondArrays() {
  bonds_.clear();
  bondsh_.clear();
}

/** Check if given angle parm exists in given angle parm array. Add if not.
  * \return Index in angle parm array.
  */
int Topology::addAngleParm(AngleParmArray& aparray, AngleParmType const& APin)
{
  // See if the AngleParm exists.
  int pidx = -1;
  for (AngleParmArray::const_iterator ap = aparray.begin();
                                      ap != aparray.end(); ++ap)
  {
    if (APin == *ap) {
      pidx = (int)(ap - aparray.begin());
      break;
    }
  }
  if (pidx == -1) {
    pidx = (int)aparray.size();
    aparray.push_back( APin );
  }
  return pidx;
}

// Topology::AddAngle() 
void Topology::AddAngle(int atom1, int atom2, int atom3, AngleParmType const& APin) {
  // See if the AngleParm exists.
  int pidx = addAngleParm( angleparm_, APin );
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
/** Check if given dihedral parm exists in given dihedral parm array. Add if not.
  * \return Index in dihedral parm array.
  */
int Topology::addTorsionParm(DihedralParmArray& dparray, DihedralParmType const& DPin)
{
  // See if the DihedralParm exists.
  int pidx = -1;
  for (DihedralParmArray::const_iterator dp = dparray.begin();
                                         dp != dparray.end(); ++dp)
  {
    if (DPin == *dp) {
      pidx = (int)(dp - dparray.begin());
      break;
    }
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
  int pidx = addTorsionParm(dihedralparm_, DPin);
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

/** \return true if any CHARMM parameters are set (based on indices). */
bool Topology::HasChamber() const {
  if (!ub_.empty()) return true;
  if (!impropers_.empty()) return true;
  if (!nonbond_.LJ14().empty()) return true;
  return false;
}

/** Add given Charmm improper with given improper parm to Charmm improper array. */
void Topology::AddCharmmImproper(DihedralType const& imp, DihedralParmType const& IPin)
{
  int pidx = addTorsionParm(improperparm_, IPin);
  if (CheckTorsionRange(imp, "CHARMM improper")) return;
  AddCharmmImproper(imp, pidx);
}

/** Add given Charmm improper with given improper parm index. */
void Topology::AddCharmmImproper(DihedralType const& impIn, int pidxIn)
{
  if (CheckTorsionRange(impIn, "CHARMM improper")) return;
  DihedralType imp = SetTorsionParmIndex(impIn, improperparm_, pidxIn, "CHARMM improper");
  // Update Charmm improper array.
  impropers_.push_back( imp );
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

/** Put all atoms in a single molecule. Mostly intended for cases
  * where you want a pseudo-topology and do not really care about
  * molecule info.
  */
int Topology::SetSingleMolecule() {
  molecules_.clear();
  molecules_.push_back( Molecule(0, Natom()) );
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

/** \return True if total mass of selected atoms is zero. */
bool Topology::MaskHasZeroMass(AtomMask const& mask) const {
  double totalMass = 0;
  for (AtomMask::const_iterator it = mask.begin(); it != mask.end(); ++it)
    totalMass += atoms_[*it].Mass();
  if (totalMass == 0) {
    mprintf("Warning: The total mass of atoms in mask '%s' is zero; cannot use mass-weighting.\n",
              mask.MaskString());
    return true;
  }
  return false;
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
        newparm.SetPk( newparm.Pk() * scale_factor ); //newparm.Pk() *= scale_factor;
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
      dk->SetPk( dk->Pk() * scale_factor ); //dk->Pk() *= scale_factor;
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
  // PDB info
  TopVecStrip<int> stripInt;
  stripInt.Strip(newParm->pdbSerialNum_, pdbSerialNum_, MapIn);
  TopVecStrip<char> stripChar;
  stripChar.Strip(newParm->atom_altloc_, atom_altloc_, MapIn);
  TopVecStrip<float> stripFloat;
  stripFloat.Strip(newParm->occupancy_, occupancy_, MapIn);
  stripFloat.Strip(newParm->bfactor_, bfactor_, MapIn);
  newParm->missingRes_ = missingRes_;
  newParm->missingHet_ = missingHet_;
 
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
    if (!nonbond_.LJ14().empty())
      //newParm->nonbond_.SetNLJ14terms( (oldTypeArray.size()*(oldTypeArray.size()+1))/2 );
      newParm->nonbond_.SetNLJ14terms( newParm->nonbond_.NBarray().size() );
    if (!nonbond_.LJC_Array().empty())
      newParm->nonbond_.SetNLJCterms( newParm->nonbond_.NBarray().size() );
    for (int a1idx = 0; a1idx != (int)oldTypeArray.size(); a1idx++)
    {
      int atm1 = oldTypeArray[a1idx];
      for (int a2idx = a1idx; a2idx != (int)oldTypeArray.size(); a2idx++)
      {
        int atm2 = oldTypeArray[a2idx];
        int oldnbidx = nonbond_.GetLJindex( atm1, atm2 );
        int newnbidx;
        if (oldnbidx > -1) {
          // This is a traditional LJ 6-12 term. Because of the way the LJ 1-4
          // code is laid out in sander/pmemd the LJ matrix has to be laid out
          // indepdendent of the nonbond index array.
          newnbidx = newParm->nonbond_.AddLJterm( a1idx, a2idx, nonbond_.NBarray(oldnbidx) );
        } else {
          // This is an old LJ 10-12 hbond term. Add one to the LJ 6-12 matrix
          // and one to the hbond since that seems to be the convention.
          newnbidx = newParm->nonbond_.AddLJterm( a1idx, a2idx, NonbondType() );
          newParm->nonbond_.AddHBterm( a1idx, a2idx, nonbond_.HBarray((-oldnbidx)-1) );
        }
        //int newnbidx = newParm->nonbond_.GetLJindex( a1idx, a2idx );
        //mprintf("DEBUG: oldtypei=%i oldtypej=%i Old NB index=%i, newtypi=%i newtypej=%i new NB idx=%i testidx=%i\n", 
        //        atm1, atm2, oldnbidx, a1idx, a2idx, newnbidx, testidx);
        if (!nonbond_.LJ14().empty()) {
          // Update LJ 1-4 as well. No need to worry about hbond terms here,
          // just recalculate the old index and determine new one.
          //int ibig = std::max(atm1, atm2) + 1;
          //int isml = std::min(atm1, atm2) + 1;
          //    oldnbidx = (ibig*(ibig-1)/2+isml)-1;
          //    ibig = a2idx + 1;
          //    isml = a1idx + 1;
          //int newnbidx = (ibig*(ibig-1)/2+isml)-1;
          newParm->nonbond_.SetLJ14( newnbidx ) = nonbond_.LJ14()[oldnbidx];
        }
        if (!nonbond_.LJC_Array().empty()) {
          // Update LJC
          newParm->nonbond_.SetLJC( newnbidx, nonbond_.LJC_Array( oldnbidx ) );
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
  if (!ff_desc_.empty())
    newParm->ff_desc_ = ff_desc_;
  if (!ub_.empty()) {
    // Urey-Bradley
    newParm->ub_ = StripBondArray( ub_, atomMap );
    parmMap.assign( ubparm_.size(), -1 ); // Map[oldidx] = newidx
    StripBondParmArray( newParm->ub_, parmMap,
                        newParm->ubparm_, ubparm_ );
  }
  if (!impropers_.empty()) {
    // Impropers
    newParm->impropers_ = StripDihedralArray(impropers_, atomMap);
    parmMap.assign( improperparm_.size(), -1 );
    StripDihedralParmArray( newParm->impropers_, parmMap,
                            newParm->improperparm_, improperparm_ );
    // NOTE 1-4 LJ parameters handled above
  }
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
  // Amber extra info.
  TopVecStrip<NameType> stripNameType;
  stripNameType.Strip(newParm->tree_, tree_, MapIn);
  stripInt.Strip(newParm->ijoin_, ijoin_, MapIn);
  stripInt.Strip(newParm->irotat_, irotat_, MapIn);
 
  // Determine number of extra points
  newParm->DetermineNumExtraPoints();

  return newParm;
}

/** Split atoms selected in a single residue into a new residue. */
int Topology::SplitResidue(AtomMask const& maskIn, NameType const& newName,
                           std::vector<int>& atomMap)
{
  if (maskIn.Nselected() == 0) {
    mprinterr("Error: SplitResidue: No atoms selected.\n");
    return 1;
  }
  int tgtResNum = Atoms()[maskIn[0]].ResNum();
  Residue const& res = residues_[tgtResNum];
  // Check that all atoms are in the same residue.
  if (maskIn.Nselected() > 1) {
    //int lastAtom = maskIn[0];
    for (int idx = 1; idx < maskIn.Nselected(); idx++) {
      //if (maskIn[idx] - lastAtom > 1) {
      //  mprinterr("Error: SplitResidue: Atoms '%s' and '%s' are not consecutive.\n",
      //            AtomMaskName(maskIn[idx]).c_str(), AtomMaskName(lastAtom).c_str());
      //  return 1;
      //}
      //lastAtom = maskIn[idx];
      if (Atoms()[maskIn[idx]].ResNum() != tgtResNum) {
        mprinterr("Error: SplitResidue: Atoms '%s' and '%s' are in different residues.\n",
                  AtomMaskName(maskIn[idx]).c_str(), AtomMaskName(maskIn[0]).c_str());
        return 1;
      }
    }
  }
  // Need to re-order the topology so that selected atoms now come at the
  // end of the residue they are a part of.
  atomMap.clear();
  atomMap.reserve(Natom());
  int r0firstAtom = -1;
  int r0lastAtom = -1;
  int r1firstAtom = -1;
  int r1lastAtom = -1;
  int newAt = 0;
  for (int at = 0; at < residues_[tgtResNum].FirstAtom(); at++, newAt++)
    atomMap.push_back(at);
  // Add unselected atoms of residue first
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
    if (!maskIn.IsSelected(at)) {
      if (r0firstAtom == -1)
        r0firstAtom = newAt;
      atomMap.push_back(at);
      newAt++;
    }
  }
  r0lastAtom = newAt;
  // Now add selected atoms
  for (AtomMask::const_iterator it = maskIn.begin(); it != maskIn.end(); ++it) {
    if (r1firstAtom == -1)
      r1firstAtom = newAt;
    atomMap.push_back(*it);
    newAt++;
  }
  r1lastAtom = newAt;
  // Add remaining atoms
  if (tgtResNum+1 < Nres()) {
    for (int at = residues_[tgtResNum+1].FirstAtom(); at != Natom(); at++)
      atomMap.push_back(at);
  }
  //mprintf("DEBUG: New res0: index %i to %i\n", r0firstAtom, r0lastAtom);
  //mprintf("DEBUG: New res1: index %i to %i\n", r1firstAtom, r1lastAtom);

  // Reorder topology
  Topology* newTop = ModifyByMap( atomMap, false );
  if (newTop == 0) {
    mprinterr("Internal Error: SplitResidue: Could not reorder the topology.\n");
    return 1;
  }
  //mprintf("DEBUG: New res0: Atoms %s to %s\n",
  //        newTop->AtomMaskName(r0firstAtom).c_str(), newTop->AtomMaskName(r0lastAtom-1).c_str());
  //mprintf("DEBUG: New res1: Atoms %s to %s\n",
  //        newTop->AtomMaskName(r1firstAtom).c_str(), newTop->AtomMaskName(r1lastAtom-1).c_str());
  //DEBUG
//  *this = *newTop; // DEBUG
//  delete newTop; // DEBUG
//  return 0; // DEBUG
  // END DEBUG
  // Decide insertion codes.
  char icode0, icode1;
  bool recode = false;
  if (res.Icode() == ' ') {
    icode0 = 'A';
    icode1 = 'B';
  } else {
    icode0 = res.Icode();
    icode1 = icode0 + 1;
    recode = true;
  }
  // Now redo the residue information
  newTop->residues_.clear();
  newTop->residues_.reserve(Nres()+1);
  // First add residues up to this one
  for (int rnum = 0; rnum < tgtResNum; rnum++)
    newTop->residues_.push_back( residues_[rnum] );
  // Add non-selected part of residue
  //mprintf("DEBUG: Last residue: %s %i - %i\n", *(newTop->residues_.back().Name()), newTop->residues_.back().FirstAtom()+1, newTop->residues_.back().LastAtom());
  newTop->residues_.push_back( Residue(residues_[tgtResNum], r0firstAtom, r0lastAtom) );
  newTop->residues_.back().SetIcode( icode0 );
  //mprintf("DEBUG: New R0: %s %i - %i\n", *(newTop->residues_.back().Name()), newTop->residues_.back().FirstAtom()+1, newTop->residues_.back().LastAtom());
  // Update atoms in selected part of residue
  for (int at = r1firstAtom; at != r1lastAtom; at++) {
    //mprintf("DEBUG: Set atom %i\n", at+1);
    newTop->atoms_[at].SetResNum( newTop->residues_.size() );
  }
  // Add selected part of residue
  newTop->residues_.push_back( Residue(residues_[tgtResNum], r1firstAtom, r1lastAtom) );
  newTop->residues_.back().SetIcode( icode1 );
  newTop->residues_.back().SetName( newName );
  //mprintf("DEBUG: New R1: %s %i - %i\n", *(newTop->residues_.back().Name()), newTop->residues_.back().FirstAtom()+1, newTop->residues_.back().LastAtom());
  int newResNum = (int)newTop->residues_.size();
  // Add remaining residues
  if (recode) {
    for (int rnum = tgtResNum+1; rnum < Nres(); rnum++, newResNum++) {
      newTop->residues_.push_back( residues_[rnum] );
      if (newTop->residues_.back().OriginalResNum() == res.OriginalResNum() &&
          newTop->residues_.back().ChainID() == res.ChainID())
        newTop->residues_.back().SetIcode(++icode1);
      for (int at = newTop->residues_.back().FirstAtom();
               at != newTop->residues_.back().LastAtom(); ++at)
        newTop->atoms_[at].SetResNum( newResNum );
    }
  } else {
    for (int rnum = tgtResNum+1; rnum < Nres(); rnum++, newResNum++) {
      newTop->residues_.push_back( residues_[rnum] );
      for (int at = newTop->residues_.back().FirstAtom();
               at != newTop->residues_.back().LastAtom(); ++at)
        newTop->atoms_[at].SetResNum( newResNum );
      //mprintf("DEBUG: Res: %s %i - %i\n", *(newTop->residues_.back().Name()), newTop->residues_.back().FirstAtom()+1, newTop->residues_.back().LastAtom());
    }
  }

  *this = *newTop;
  delete newTop;

  return 0;
}

/** Split atoms selected in a single residue into a new residue. */
int Topology::SplitResidue(AtomMask const& maskIn, NameType const& newName)
{
  std::vector<int> atomMap;
  return SplitResidue(maskIn, newName, atomMap);
}

// -----------------------------------------------------------------------------
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
            //if (olddih->Type() != DihedralType::NORMAL && (newA3 == 0 || newA4 == 0))
            //  dihOut.push_back( DihedralType( newA4, newA3, newA2, newA1, 
            //                                  olddih->Type(), olddih->Idx() ) );
            //else
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

// -----------------------------------------------------------------------------
/** Loop over atom type pairs, determine what the A and B
  * parameters should be based on the LJ radius and depth,
  * then compare to the existing A and B parameters.
  * Warn if it does not match, store the current A and B
  * parameters.
  * \param atomTypesOut The atom types containing radii and depths.
  * \param LJ612out The output LJ 6-12 A and B paramters.
  * \param LJ1012ptr If not 0, the output LJ 10-12 parameters.
  */
/*static inline void GetLJterms(ParmHolder<AtomType> const& atomTypesOut,
                              ParmHolder<NonbondType>& LJ612out,
                              ParmHolder<HB_ParmType>* LJ1012ptr,
                              NonbondParmType const& NB0,
                              bool isLJ14,
                              int debugIn)
{
    // Do atom type pairs, check for off-diagonal elements.
    // Explicitly store pairs instead of regenerating to avoid round-off issues.
    unsigned int nModifiedOffDiagonal = 0;
    for (ParmHolder<AtomType>::const_iterator i1 = atomTypesOut.begin(); i1 != atomTypesOut.end(); ++i1)
    {
      for (ParmHolder<AtomType>::const_iterator i2 = i1; i2 != atomTypesOut.end(); ++i2)
      {
        NameType const& name1 = i1->first[0];
        NameType const& name2 = i2->first[0];
        TypeNameHolder types(2);
        types.AddName( name1 );
        types.AddName( name2 );
        // Extract original nonbonded parameters for this type pair.
        AtomType const& type1 = i1->second;
        AtomType const& type2 = i2->second;
        int idx1 = type1.OriginalIdx();
        int idx2 = type2.OriginalIdx();
        int idx = NB0.GetLJindex( idx1, idx2 );
        if (idx < 0) {
          // This is LJ 10-12.
          //mprinterr("Error: No off-diagonal LJ for  %s %s (%i %i)\n",
          //          *name1, *name2, idx1, idx2);
          //return;
          if (LJ1012ptr != 0) {
            mprintf("DEBUG: LJ 10-12 parameters detected for %s %s (%i %i)\n",
                    *name1, *name2, idx1, idx2);
            LJ1012ptr->AddParm( types, NB0.HBarray((-idx)-1), false );
          }
        } else {
          // This is LJ 6-12.
          // Determine what A and B parameters would be.
          NonbondType lj0 = type1.LJ().Combine_LB( type2.LJ() );

          NonbondType lj1;
          if (isLJ14)
            lj1 = NB0.LJ14( idx );
          else
            lj1 = NB0.NBarray( idx );
          // Compare them
          if (lj0 != lj1) {
            nModifiedOffDiagonal++;
            //if (debugIn > 0) {
              double deltaA = fabs(lj0.A() - lj1.A());
              double deltaB = fabs(lj0.B() - lj1.B());
              mprintf("DEBUG: Potential off-diagonal LJ: %s %s expect A=%g B=%g, actual A=%g B=%g\n",
                      *name1, *name2, lj0.A(), lj0.B(), lj1.A(), lj1.B());
              mprintf("DEBUG:\tdeltaA= %g    deltaB= %g\n", deltaA, deltaB);
              double pe_a = (fabs(lj0.A() - lj1.A()) / lj0.A());
              double pe_b = (fabs(lj0.B() - lj1.B()) / lj0.B());
              mprintf("DEBUG:\tPEA= %g  PEB= %g\n", pe_a, pe_b);
            //}
          }
          LJ612out.AddParm( types, lj1, false );
        }
      } // END inner loop over atom types
    } // END outer loop over atom types
    if (nModifiedOffDiagonal > 0)
      mprintf("Warning: %u modified off-diagonal LJ terms present.\n", nModifiedOffDiagonal);
}*/



// -----------------------------------------------------------------------------
/** Append a topology to another with default options */
int Topology::AppendTop(Topology const& NewTop) {
  return AppendTop(NewTop, 0, false, false);
}

/** Append a topology to another */
int Topology::AppendTop(Topology const& NewTop, int verbose_, bool reduce_bond_params, bool reduce_angle_params) {
  using namespace Cpptraj::Parm;
  Topology& topOut = *this;

  unsigned int atomOffset = (unsigned int)topOut.Natom();
  unsigned int molOffset = (unsigned int)topOut.Nmol();
  if (debug_ > 0)
    mprintf("DEBUG: Appending '%s' to '%s' (atom offset= %u mol offset= %u)\n",
            NewTop.c_str(), topOut.c_str(), atomOffset, molOffset);
  //int resOffset = (int)residues_.size();

  // Save nonbonded parameters from each topology
  ParmHolder<AtomType> myAtomTypes, newAtomTypes;
  ParmHolder<NonbondType> myNB, newNB;
  ParmHolder<NonbondType> my14, new14;
  ParmHolder<HB_ParmType> myHB, newHB;
  ParmHolder<double> myLJC, newLJC;
  Cpptraj::Parm::GetParams GP;
  GP.SetDebug( debug_ );
  GP.GetLJAtomTypes( myAtomTypes, myNB, my14, myLJC, myHB, topOut.Atoms(), topOut.Nonbond() );
  GP.GetLJAtomTypes( newAtomTypes, newNB, new14, newLJC, newHB, NewTop.Atoms(), NewTop.Nonbond() );
  // Create combined nonbond parameter set
  int nAtomTypeUpdated = UpdateParameters< ParmHolder<AtomType> >( myAtomTypes, newAtomTypes, "atom type", verbose_ );
  int nLJparamsUpdated = UpdateParameters< ParmHolder<NonbondType> >( myNB, newNB, "LJ 6-12", verbose_ );
  int n14paramsUpdated = UpdateParameters< ParmHolder<NonbondType> >( my14, new14, "LJ 6-12 1-4", verbose_ );
  int nljcparamsUpdated = UpdateParameters< ParmHolder<double> >( myLJC, newLJC, "LJ C", verbose_ );
  int nHBparamsUpdated = UpdateParameters< ParmHolder<HB_ParmType> >( myHB, newHB, "LJ 10-12", verbose_ );
  if (verbose_ > 0)
    mprintf("\t%i atom types updated, %i LJ 6-12 params updated, %i 1-4 params updated, %i LJC params updated, %i LJ 10-12 params updated.\n",
            nAtomTypeUpdated, nLJparamsUpdated, n14paramsUpdated, nljcparamsUpdated, nHBparamsUpdated);

  // Add incoming topology bond/angle/dihedral/cmap arrays to this one.
  MergeBondArrays(reduce_bond_params, topOut.ModifyBonds(), topOut.ModifyBondsH(), topOut.ModifyBondParm(), topOut.Atoms(),
                  NewTop.Bonds(), NewTop.BondsH(), NewTop.BondParm(), NewTop.Atoms());
  MergeAngleArrays(reduce_angle_params, topOut.ModifyAngles(), topOut.ModifyAnglesH(), topOut.ModifyAngleParm(), topOut.Atoms(),
                   NewTop.Angles(), NewTop.AnglesH(), NewTop.AngleParm(), NewTop.Atoms());
  MergeDihedralArrays(topOut.ModifyDihedrals(), topOut.ModifyDihedralsH(), topOut.ModifyDihedralParm(), topOut.Atoms(),
                      NewTop.Dihedrals(), NewTop.DihedralsH(), NewTop.DihedralParm(), NewTop.Atoms());
  MergeCmapArrays(topOut.ModifyCmap(), topOut.ModifyCmapGrid(), topOut.Atoms(), topOut.Residues(),
                  NewTop.Cmap(), NewTop.CmapGrid(), NewTop.Atoms(), NewTop.Residues());
  MergeBondArray(topOut.ModifyUB(), topOut.ModifyUBparm(), topOut.Atoms(),
                 NewTop.UB(), NewTop.UBparm(), NewTop.Atoms());
  MergeImproperArray(topOut.ModifyImpropers(), topOut.ModifyImproperParm(), topOut.Atoms(),
                     NewTop.Impropers(), NewTop.ImproperParm(), NewTop.Atoms());

  // Append incoming atoms to this topology.
  for (AtArray::const_iterator atom = NewTop.begin(); atom != NewTop.end(); ++atom)
  {
    if (debug_ > 1)
      mprintf("DEBUG: %6li %s %s %4i\n", atom-NewTop.begin(),
              *(atom->Name()), *(atom->Type()), atom->TypeIndex());
    Atom CurrentAtom = *atom;
    Residue const& res = NewTop.Res( CurrentAtom.ResNum() );
    // Bonds need to be cleared and re-added.
    CurrentAtom.ClearBonds();
    for (Atom::bond_iterator bat = atom->bondbegin(); bat != atom->bondend(); ++bat)
      CurrentAtom.AddBondToIdx( *bat + atomOffset );

    topOut.addTopAtom( CurrentAtom,
                       Residue(res.Name(), res.OriginalResNum(), res.Icode(), res.ChainID()),
                       atom->MolNum()+molOffset, NewTop.Mol(atom->MolNum()).IsSolvent() );
  }

  // EXTRA ATOM INFO
  TopVecAppend<NameType> appendNameType;
  appendNameType.Append( topOut.ModifyTreeChainClassification(), NewTop.TreeChainClassification(), NewTop.Natom() );
  TopVecAppend<int> appendInt;
  appendInt.Append( topOut.ModifyJoinArray(), NewTop.JoinArray(), NewTop.Natom() );
  appendInt.Append( topOut.ModifyRotateArray(), NewTop.RotateArray(), NewTop.Natom() );
  appendInt.Append( topOut.ModifyPdbSerialNum(), NewTop.PdbSerialNum(), NewTop.Natom() );
  TopVecAppend<char> appendChar;
  appendChar.Append( topOut.ModifyAtomAltLoc(), NewTop.AtomAltLoc(), NewTop.Natom() );
  TopVecAppend<float> appendFloat;
  appendFloat.Append( topOut.ModifyOccupancy(), NewTop.Occupancy(), NewTop.Natom() );
  appendFloat.Append( topOut.ModifyBfactor(), NewTop.Bfactor(), NewTop.Natom() );

  // Need to regenerate nonbonded info
  if (verbose_ > 0)
    mprintf("\tRegenerating nonbond parameters.\n");
  AssignParams assign;
  assign.SetDebug( debug_ );
  assign.SetVerbose( verbose_ );
  assign.AssignNonbondParams( topOut, myAtomTypes, myNB, my14, myLJC, myHB );

  // The version of AddTopAtom() with molecule number already determines
  // molecules and number of solvent molecules.
  // Just need to determine the number of extra points.
  topOut.DetermineNumExtraPoints();

  // GB radii set string
  if (!NewTop.GBradiiSet().empty()) {
    if (topOut.GBradiiSet().empty())
      topOut.SetGBradiiSet( NewTop.GBradiiSet() );
    else {
      // Do not repeat a GB radius string
      std::size_t pos = topOut.GBradiiSet().find( NewTop.GBradiiSet() );
      if (pos == std::string::npos) {
        std::string newName = topOut.GBradiiSet() + "+" + NewTop.GBradiiSet();
        if (newName.size() > 80) {
          mprintf("Warning: New radius set name is > 80 characters: '%s'\n", newName.c_str());
          mprintf("Warning: This will be truncated to 80 characters in an Amber Topology.\n");
        }
        topOut.SetGBradiiSet( newName );
      }
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** \return True if any atom has a non-zero charge. */
bool Topology::HasChargeInfo() const {
  for (std::vector<Atom>::const_iterator at = atoms_.begin();
                                         at != atoms_.end(); ++at)
    if (at->Charge() > 0.0 || at->Charge() < 0.0)
      return true;
  return false;
}

/** \return Total charge */
double Topology::TotalCharge() const {
  double sumQ = 0.0;
  for (std::vector<Atom>::const_iterator at = atoms_.begin();
                                         at != atoms_.end(); ++at)
    sumQ += at->Charge();
  return sumQ;
}

/** Redistribute charge on atoms to match given total target charge. */
int Topology::RedistributeCharge(double charge) {
  //mprintf("DEBUG: Redistribute charge for %s, total charge = %g\n", topIn.c_str(), charge);
  double pcharge = 0;
  double ncharge = 0;
  for (unsigned int iat = 0; iat != atoms_.size(); iat++) {
    if (atoms_[iat].Charge() > 0)
      pcharge += atoms_[iat].Charge();
    else if (atoms_[iat].Charge() < 0)
      ncharge += atoms_[iat].Charge();
  }
  //if (fabs(pcharge) < Constants::SMALL)
  bool PchargeZero = false;
  if (pcharge == 0) {
    mprintf("\tTotal positive charge is 0.0\n");
    PchargeZero = true;
  }
  bool NchargeZero = false;
  //if (fabs(ncharge) < Constants::SMALL)
  if (ncharge == 0) {
    mprintf("\tTotal negative charge is 0.0\n");
    NchargeZero = true;
  }
  if (!PchargeZero && !NchargeZero) {
    //double total_charge = 0;
    for (unsigned int iat = 0; iat != atoms_.size(); iat++) {
      double delta = atoms_[iat].Charge() * (charge - pcharge - ncharge) / (pcharge - ncharge);
      if (atoms_[iat].Charge() >= 0) {
        atoms_[iat].SetCharge( atoms_[iat].Charge() + delta );
      } else {
        atoms_[iat].SetCharge( atoms_[iat].Charge() - delta );
      }
      //total_charge += topIn[iat].Charge();
    }
    //mprintf("DEBUG: Total charge after redistribute: %g\n", total_charge);
  }
  return 0;
}
