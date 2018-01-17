#include "AtomTypeArray.h"
#include "CpptrajStdio.h"

/** \return true if type name already present. */
bool AtomTypeArray::AddAtomType(NameType const& nameIn, AtomType const& typeIn)
{
  // See if already present.
  Tmap::iterator it = nameToIdx_.find( nameIn );
  if (it == nameToIdx_.end()) {
    // New type
    int idx = types_.size();
    types_.push_back( typeIn );
    nameToIdx_.insert( Tpair(nameIn, idx) );
    mprintf("\tAdded atom type '%s'\n", *nameIn);
    return false;
  }
  mprintf("\tType '%s' already present.\n");
  return true;
}
