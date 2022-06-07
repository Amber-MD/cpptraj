#include "StructureRoutines.h"
#include "../Atom.h"
#include "../CpptrajStdio.h"
#include "../Residue.h"

/** Global Structure debug level. */
static int structure_debug_ = 0;

using namespace Cpptraj::Structure;

/// Set the Structure debug level
void SetStructureDebugLevel(int debugIn) { structure_debug_ = debugIn; }

/// \return the Structure debug level
int StructureDebugLevel() { return structure_debug_; }

/// Used to change residue name to nameIn
void ChangeResName(Residue& res, NameType const& nameIn) {
  if (res.Name() != nameIn) {
    if (structure_debug_ > 0) mprintf("\t    Changing residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

/// Used to change atom name to nameIn
void ChangeAtomName(Atom& atm, NameType const& nameIn) {
  if (atm.Name() != nameIn) {
    if (structure_debug_ > 0) mprintf("\t    Changing atom %s to %s\n", *(atm.Name()), *nameIn);
    atm.SetName( nameIn );
  }
}


