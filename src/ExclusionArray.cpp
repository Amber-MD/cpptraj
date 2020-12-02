#include "ExclusionArray.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

/** CONSTRUCTOR */
ExclusionArray::ExclusionArray() {}

/** Set up exclusion array from given list of atoms and atom mask. */
int ExclusionArray::SetupExcluded(std::vector<Atom> const& atoms, AtomMask const& maskIn)
{
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
    // Always exclude self
    Excluded_[idx].insert( idx );
    int at = maskIn[idx];
    for (Atom::excluded_iterator excluded_atom = atoms[at].excludedbegin();
                                 excluded_atom != atoms[at].excludedend();
                               ++excluded_atom)
    {
      if (Cmask.AtomInCharMask(*excluded_atom))
      {
        // Find excluded atoms index in maskIn
        int excluded_idx = atToIdx[*excluded_atom];
        Excluded_[idx         ].insert( excluded_idx );
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
