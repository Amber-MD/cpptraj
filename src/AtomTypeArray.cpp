#include "AtomTypeArray.h"
#include "CpptrajStdio.h"

/** Add atom type to array if not already present.
  * \return true if type name already present.
  */
bool AtomTypeArray::AddAtomType(NameType const& nameIn, AtomType const& typeIn)
{
  // See if already present.
  Tmap::iterator it = nameToIdx_.find( nameIn );
  if (it == nameToIdx_.end()) {
    // New type
    int idx = types_.size();
    types_.push_back( typeIn );
    nameToIdx_.insert( Tpair(nameIn, idx) );
    if (debug_ > 0)
      mprintf("\tAdded atom type '%s', mass=%f radius=%f depth=%f\n", *nameIn,
              typeIn.Mass(), typeIn.LJ().Radius(), typeIn.LJ().Depth());
    return false;
  }
  if (debug_ > 0) mprintf("\tType '%s' already present.\n", *nameIn);
  return true;
}

/** Add new or find existing atom type based on name.
  * \return Index of new/added atom type.
  */
int AtomTypeArray::CheckForAtomType(NameType const& nameIn, AtomType const& typeIn) {
  int idx = -1;
  Tmap::iterator it = nameToIdx_.find( nameIn );
  if (it == nameToIdx_.end()) {
    // New placholder type
    idx = types_.size();
    types_.push_back( typeIn );
    nameToIdx_.insert( Tpair(nameIn, idx) );
    if (debug_ > 0) mprintf("\tAdded placholder atom type '%s' (%i)\n", *nameIn, idx);
  } else {
    idx = it->second;
    if (debug_ > 0) mprintf("\tType '%s' already present (%i).\n", *nameIn, idx);
  }
  return idx;
}

// AtomTypeArray::CheckForAtomType()
int AtomTypeArray::CheckForAtomType(NameType const& nameIn) {
  return CheckForAtomType( nameIn, AtomType() );
}

/** \return Atom type index if given atom type name is present, -1 otherwise. */
int AtomTypeArray::AtomTypeIndex(NameType const& nameIn) {
  int idx = -1;
  Tmap::iterator it = nameToIdx_.find( nameIn );
  if (it != nameToIdx_.end())
    idx = it->second;
  return idx;
}
