#include <algorithm> // sort
#include "AtomMap.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AtomMap::AtomMap() :
  debug_(0)
{}

MapAtom& AtomMap::operator[](int idx) {
  if (idx < 0 || idx >= (int)mapatoms_.size()) {
    mprinterr("Error: AtomMap::operator[]: Index %i out of range.\n",idx);
    return EMPTYMAPATOM;
  }
  return mapatoms_[idx];
}

int AtomMap::Natom() {
  return (int)mapatoms_.size();
}

void AtomMap::SetDebug(int debugIn) {
  debug_ = debugIn;
}

int AtomMap::Setup(Topology *TopIn) {
  // Copy atoms
  mapatoms_.clear();
  for (Topology::atom_iterator atom = TopIn->begin(); atom != TopIn->end(); atom++) {
    // This sets up 1 char atom name based on atom element
    mapatoms_.push_back( *atom );
  }
  return 0;
}

void AtomMap::ResetMapping() {
  for (std::vector<MapAtom>::iterator matom = mapatoms_.begin();
                                      matom != mapatoms_.end(); matom++)
  {
    (*matom).SetNotMapped();
    (*matom).SetNotComplete();
  }
}

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

void AtomMap::MarkAtomComplete(int atom) {
  if (atom<0 || atom >= (int)mapatoms_.size()) return;
  if (!mapatoms_[atom].IsMapped() && debug_==0) return;
  if ( mapatoms_[atom].Complete() && debug_==0) return;
  int nunique = 0;
  for (Atom::bond_iterator bondedAtom = mapatoms_[atom].bondbegin();
                           bondedAtom != mapatoms_[atom].bondend(); bondedAtom++)
    if (mapatoms_[*bondedAtom].IsMapped())
      ++nunique;
  if (mapatoms_[atom].IsUnique() && nunique==mapatoms_[atom].Nbonds())
    mapatoms_[atom].SetComplete();
  if (debug_>0) {
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

void AtomMap::CheckForCompleteAtoms() {
  for (int atom = 0; atom < (int)mapatoms_.size(); atom++)
    MarkAtomComplete(atom);
}

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
      if (N_single_atoms<3)
        (*matom).SetChiral();
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
      if ((*matom).IsChiral()) mprintf(" CHIRAL!");
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
