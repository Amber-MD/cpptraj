#include <algorithm> // sort
#include "AtomMap.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AtomMap::AtomMap() :
  debug_(0)
{}

/// Blank AtomMap for empty return of bracket operator
MapAtom AtomMap::EMPTYMAPATOM = MapAtom();

// AtomMap::operator[]()
/** Return a reference to the specifed MapAtom */
MapAtom& AtomMap::operator[](int idx) {
  if (idx < 0 || idx >= (int)mapatoms_.size()) {
    mprinterr("Error: AtomMap::operator[]: Index %i out of range.\n",idx);
    return EMPTYMAPATOM;
  }
  return mapatoms_[idx];
}

// AtomMap::Natom()
/// Return the number of atoms in the AtomMap.
int AtomMap::Natom() {
  return (int)mapatoms_.size();
}

// AtomMap::SetDebug()
/// Set the debug level of the AtomMap.
void AtomMap::SetDebug(int debugIn) {
  debug_ = debugIn;
}

/// Check if 1 char name set to 0, means unidentified element.
bool AtomMap::InvalidElement() {
  if (mapatoms_.back().CharName() == 0) {
    mprinterr("Error: AtomMap: Mapping currently not supported for element %s\n",
              mapatoms_.back().ElementName());
    return true;
  }
  return false;
}

// AtomMap::Setup()
/** Copy all atoms from input topology to this AtomMap. */
int AtomMap::Setup(Topology const& TopIn) {
  mapatoms_.clear();
  for (Topology::atom_iterator atom = TopIn.begin(); atom != TopIn.end(); atom++) {
    // This sets up 1 char atom name based on atom element
    mapatoms_.push_back( *atom );
    if (InvalidElement()) return 1;
  }
  return 0;
}

int AtomMap::SetupResidue(Topology const& topIn, int resnum) {
  mapatoms_.clear();
  int firstAtom = topIn.Res(resnum).FirstAtom();
  int lastAtom = topIn.Res(resnum).LastAtom();
  //mprintf("DEBUG:\tResidue %i, atoms %i to %i\n", resnum + 1, firstAtom+1, lastAtom);
  for (int atom = firstAtom; atom < lastAtom; ++atom) {
    mapatoms_.push_back( topIn[atom] );
    if (InvalidElement()) return 1;
    // Add bonds for this residue 
    mapatoms_.back().ClearBonds();
    for (Atom::bond_iterator bndatm = topIn[atom].bondbegin();
                             bndatm != topIn[atom].bondend(); ++bndatm)
    {
      //mprintf("DEBUG:\t\tOriginal bond %u-%i", atom+1, *bndatm+1);
      if (*bndatm >= firstAtom && *bndatm < lastAtom) { 
        int newbndatm = *bndatm - firstAtom;
        mapatoms_.back().AddBond(newbndatm);
        //mprintf(", new bond %i-%i", mapatoms_.size(), newbndatm+1);
      }
      //mprintf("\n");
    }
  }
  return 0;
}

// AtomMap::ResetMapping()
/** Reset any previously set mapping information. */
void AtomMap::ResetMapping() {
  for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                      matom != mapatoms_.end(); matom++)
  {
    (*matom).SetNotMapped();
    (*matom).SetNotComplete();
  }
}

// AtomMap::BondIsRepeated()
/** Check if the atomID of the specified atom (bondedAtom) bonded to <atom> 
  * is the same as the atomID of any other non-mapped atom bonded to <atom>.
  */
bool AtomMap::BondIsRepeated(int atom, int bondedAtom) {
  // If 1 or no bonds, atom cant possibly be repeated
  if (mapatoms_[atom].Nbonds() > 1) {
    for (Atom::bond_iterator bondedAtom2 = mapatoms_[atom].bondbegin();
                             bondedAtom2 != mapatoms_[atom].bondend(); bondedAtom2++)
    {
      if (mapatoms_[*bondedAtom2].IsMapped()) continue;
      if (mapatoms_[bondedAtom].AtomID() == mapatoms_[*bondedAtom2].AtomID())
        return true;
    }
  }
  return false;
}

// AtomMap::DetermineAtomIDs()
/** Give each atom an identifier (atomID) based on what atoms are bonded to 
  * it. The first part of the atomID is the atom itself, followed by an 
  * alphabetized list of bonded atoms. So C in O=C-H2 would be CHHO.
  * Then create a 'unique string' comprised of this atomID and the 
  * atomIDs of all bonded atoms (sorted). Once that is done determine
  * which unique strings are actually unique (i.e. they are not repeated
  * in this map). 
  */
void AtomMap::DetermineAtomIDs() {
  // Determine self IDs
  if (debug_>0) mprintf("ATOM IDs:\n");
  unsigned int anum = 1;
  for (std::vector<MapAtom>::iterator matom = mapatoms_.begin(); 
                                      matom != mapatoms_.end(); matom++)
  {
    std::string atomID;
    for (Atom::bond_iterator bondedAtom = (*matom).bondbegin();
                             bondedAtom != (*matom).bondend(); bondedAtom++)
    {
      atomID += mapatoms_[ *bondedAtom ].CharName();
    }
    // Sort atom ID
    sort( atomID.begin(), atomID.end() );
    // Place current atom 1 char name at beginning
    atomID = (*matom).CharName() + atomID;
    (*matom).SetAtomID( atomID );
    if (debug_>0) mprintf("  Atom %u : %s\n",anum, atomID.c_str());
    ++anum;
  }
  
  // Create a unique ID for each atom based on Atom IDs
  for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                      matom != mapatoms_.end(); matom++)
  {
    std::string unique = (*matom).AtomID();
    for (Atom::bond_iterator bondedAtom = (*matom).bondbegin();
                             bondedAtom != (*matom).bondend(); bondedAtom++)
    {
      unique += mapatoms_[ *bondedAtom ].AtomID();
    }
    sort( unique.begin(), unique.end() );
    // NOTE: SetUnique also resets the dupliciated counter.
    (*matom).SetUnique( unique );
  }

  // Determine which unique IDs are duplicated - set isUnique flag
  for (unsigned int i = 0; i < mapatoms_.size()-1; i++) {
    for (unsigned int j = i+1; j < mapatoms_.size(); j++) {
      if ( mapatoms_[i].Unique() == mapatoms_[j].Unique() ) {
        // This unique string is duplicated
        mapatoms_[i].IsDuplicated();
        mapatoms_[j].IsDuplicated();
      }
    }
  }

  // DEBUG
  if (debug_ > 0) {
    mprintf("UNIQUE IDs:\n");
    anum = 1;
    for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                        matom != mapatoms_.end(); matom++)
    {
      mprintf("  Atom %6u [%3i]: %s",anum,(*matom).Nduplicated(),(*matom).Unique().c_str());
      if ((*matom).IsUnique()) mprintf(" UNIQUE!");
      mprintf("\n");
      ++anum;
    }
  }
}

// AtomMap::MarkAtomComplete()
/** If atom is mapped and all bonded atoms are mapped mark atom as completely 
  * mapped.
  * If printAtoms is true print isMapped value for this atom and all atoms
  * bonded to it.
  */
void AtomMap::MarkAtomComplete(int atom, bool printAtoms) {
  if (atom<0 || atom >= (int)mapatoms_.size()) return;
  if (!mapatoms_[atom].IsMapped() && !printAtoms) return;
  if ( mapatoms_[atom].Complete() && !printAtoms) return;
  int nunique = 0;
  for (Atom::bond_iterator bondedAtom = mapatoms_[atom].bondbegin();
                           bondedAtom != mapatoms_[atom].bondend(); bondedAtom++)
    if (mapatoms_[*bondedAtom].IsMapped())
      ++nunique;
  if (mapatoms_[atom].IsUnique() && nunique==mapatoms_[atom].Nbonds())
    mapatoms_[atom].SetComplete();
  if (printAtoms) {
    mprintf("  Atom %4i: %c-%1i |",atom+1,mapatoms_[atom].c_str(),
            (int)mapatoms_[atom].IsMapped());
    for (Atom::bond_iterator bondedAtom = mapatoms_[atom].bondbegin();
                           bondedAtom != mapatoms_[atom].bondend(); bondedAtom++)
    {
      mprintf(" %4i:%c-%1i",*bondedAtom+1,mapatoms_[*bondedAtom].c_str(),
              (int)mapatoms_[*bondedAtom].IsMapped());
    }
    if (mapatoms_[atom].Complete())
      mprintf(" Atom is completely mapped.");
    mprintf("\n");
  }
}

// AtomMap::CheckForCompleteAtoms()
/** Go through each atom in the map. If the atom is unique and all bonded
  * atoms are unique mark the atom as completely mapped.
  */
void AtomMap::CheckForCompleteAtoms() {
  bool printAtoms = (debug_ > 0);
  for (int atom = 0; atom < (int)mapatoms_.size(); atom++)
    MarkAtomComplete(atom,printAtoms);
}

// AtomMap::CheckBonds()
/** Checks that bonding information is present. Also checks for potential
  * chiral centers.
  */
int AtomMap::CheckBonds() {
  int total_bonds = 0;
  // Search for chiral centers by number of bonds
  for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                      matom != mapatoms_.end(); matom++)
  {
    // Sort the bonded atoms array by atom #
    (*matom).SortBonds();
    total_bonds += (*matom).Nbonds();
    if ((*matom).Nbonds() == 4) {
      // If >=3 bonds to single atoms, not chiral (e.g. -CH3)
      int N_single_atoms=0; // Count # bonds to single atoms
      for (Atom::bond_iterator bondedAtom = (*matom).bondbegin();
                               bondedAtom != (*matom).bondend(); bondedAtom++)
      {
        if (mapatoms_[*bondedAtom].Nbonds() == 1)
          ++N_single_atoms;
      }
      if (N_single_atoms<3) {
        (*matom).SetChiral();
        for (Atom::bond_iterator bondedAtom = (*matom).bondbegin();
                               bondedAtom != (*matom).bondend(); bondedAtom++)
          mapatoms_[*bondedAtom].SetBoundToChiral();
      }
    }
  }
  if (total_bonds == 0) {
    mprinterr("Error: No bond information present, required by AtomMap.\n");
    return 1;
  }

  // DEBUG
  if (debug_>0) {
    mprintf("AtomMap: Atom Bond information.\n");
    unsigned int anum = 1;
    for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                        matom != mapatoms_.end(); matom++)
    {
      mprintf("  Atom %s(%c)_%i has %i bonds.",(*matom).c_str(),(*matom).CharName(),
              anum, (*matom).Nbonds());
      if ((*matom).IsChiral()) mprintf(" CHIRAL");
      if ((*matom).BoundToChiral()) mprintf(" BOUND TO CHIRAL");
      mprintf("\n");
      for (Atom::bond_iterator bondedAtom = (*matom).bondbegin();
                               bondedAtom != (*matom).bondend(); bondedAtom++)
      {
        mprintf("    to %s(%c)_%i\n",mapatoms_[*bondedAtom].c_str(),
                mapatoms_[*bondedAtom].CharName(), *bondedAtom+1);
      }
    }
  }
  return 0;
}
