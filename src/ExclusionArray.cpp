#include "ExclusionArray.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

/** CONSTRUCTOR */
ExclusionArray::ExclusionArray() {}

/** Determine the through-bond distance from original atom to this atom; if it is less
  * than the target distance, mark it as excluded.
  */
void ExclusionArray::AtomDistance(std::vector<Atom> const& atoms,
                                  int originalAtom, int atom, int dist, ExListType& excluded,
                                  int TgtDist)
{
  // If this atom is already too far away return
  if (dist==TgtDist) return;
  // dist is less than 4 and this atom greater than original, add exclusion
  if (atom > originalAtom)
    excluded.insert( atom ); 
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = atoms[atom].bondbegin();
                           bondedatom != atoms[atom].bondend();
                           bondedatom++)
    AtomDistance(atoms, originalAtom, *bondedatom, dist+1, excluded, TgtDist);
}

/** For each atom, determine which atoms with greater atom# are within
  * TgtDist bonds (and therefore should be excluded from a non-bonded calc).
  */
/*
void ExclusionArray::DetermineExcludedAtoms(ExListType& excluded_i,
                                            std::vector<Atom> const& atoms, int atomi,
                                            int TgtDist)
{
  // A set is used since it automatically sorts itself and rejects duplicates.
  //ExListType excluded_i;
  //int natom = (int)atoms.size();
  //for (int atomi = 0; atomi < natom; atomi++) {
  //  excluded_i.clear();
    //mprintf("    Determining excluded atoms for atom %i\n",atomi+1);
    // AtomDistance recursively sets each atom bond distance from atomi
    AtomDistance(atoms, atomi, atomi, 0, excluded_i, TgtDist);
    //Excluded_.push_back( excluded_i );
    // DEBUG
    //mprintf("\tAtom %i Excluded:",atomi+1);
    //for (Atom::excluded_iterator ei = atoms_[atomi].excludedbegin(); 
    //                             ei != atoms_[atomi].excludedend(); ++ei)
    //  mprintf(" %i",*ei + 1);
    //mprintf("\n");
  //} // END loop over atomi
}*/

/** Set up exclusion array from given list of atoms and atom mask. */
int ExclusionArray::SetupExcluded(std::vector<Atom> const& atoms, AtomMask const& maskIn,
                                  int TgtDist, SelfOpt selfOpt, ListOpt listOpt)
{
  bool exclude_self;
  if (selfOpt == EXCLUDE_SELF)
    exclude_self = true;
  else // NO_EXCLUDE_SELF
    exclude_self = false;
  Excluded_.clear();
  Excluded_.resize( maskIn.Nselected() );
  // Create a character mask so we can see if atoms in excluded lists are
  // also selected.
  CharMask Cmask(maskIn.ConvertToCharMask(), maskIn.Nselected());
  // Create a map of atom number to maskIn index.
  int selectedIdx = 0;
  std::vector<int> atToIdx( Cmask.Natom(), -1 );
  for (int cidx = 0; cidx != Cmask.Natom(); cidx++)
    if (Cmask.AtomInCharMask(cidx))
      atToIdx[cidx] = selectedIdx++;
  // Loop over selected atoms
  for (int idx = 0; idx != maskIn.Nselected(); idx++)
  {
    // Exclude self if specified
    if (exclude_self) Excluded_[idx].insert( idx );
    int at = maskIn[idx];
    // Find excluded atoms for this atom.
    // Use a separate list so we can exclude via the atom mask.
    ExListType excluded_i;
    AtomDistance(atoms, at, at, 0, excluded_i, TgtDist);
    //for (Atom::excluded_iterator excluded_atom = atoms[at].excludedbegin();
    //                             excluded_atom != atoms[at].excludedend();
    //                           ++excluded_atom)
    for (ExListType::const_iterator excluded_atom = excluded_i.begin();
                                    excluded_atom != excluded_i.end();
                                  ++excluded_atom)
    {
      if (Cmask.AtomInCharMask(*excluded_atom))
      {
        // Find excluded atoms index in maskIn
        int excluded_idx = atToIdx[*excluded_atom];
        Excluded_[idx         ].insert( excluded_idx );
        if (listOpt == FULL)
          Excluded_[excluded_idx].insert( idx          );
      }
    }
  }
  unsigned int ex_size = 0;
  for (ExArrayType::const_iterator it = Excluded_.begin(); it != Excluded_.end(); ++it)
    ex_size += it->size();
  mprintf("\tMemory used by full exclusion list: %s\n",
          ByteString(ex_size * sizeof(int), BYTE_DECIMAL).c_str());
  return 0;
}

/** Set up exclusion array for all atoms. */
int ExclusionArray::SetupExcluded(std::vector<Atom> const& atoms,
                                  int TgtDist, SelfOpt selfOpt, ListOpt listOpt)
{
  AtomMask allAtoms(0, atoms.size());
  return SetupExcluded(atoms, allAtoms, TgtDist, selfOpt, listOpt);
}
