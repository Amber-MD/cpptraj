#include "StructureRoutines.h"
#include "../Atom.h"
#include "../CpptrajStdio.h"
#include "../Residue.h"

/** Global Structure debug level. */
static int structure_debug_ = 0;

/// Set the Structure debug level
void Cpptraj::Structure::SetStructureDebugLevel(int debugIn) { structure_debug_ = debugIn; }

/// \return the Structure debug level
int Cpptraj::Structure::StructureDebugLevel() { return structure_debug_; }

/// Used to change residue name to nameIn
void Cpptraj::Structure::ChangeResName(Residue& res, NameType const& nameIn) {
  if (res.Name() != nameIn) {
    if (structure_debug_ > 0) mprintf("\t    Changing residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

/// Used to change atom name to nameIn
void Cpptraj::Structure::ChangeAtomName(Atom& atm, NameType const& nameIn) {
  if (atm.Name() != nameIn) {
    if (structure_debug_ > 0) mprintf("\t    Changing atom %s to %s\n", *(atm.Name()), *nameIn);
    atm.SetName( nameIn );
  }
}


