#include "ParameterHolders.h"
#include "CpptrajStdio.h"

int BondParmHolder::AddBondParm(AtomTypeHolder const& types, BondParmType const& bp,
                                bool allowUpdate)
{
  if (types.Size() != 2) {
    mprinterr("Internal Error: BondParmHolder::AddBondParm(): # types is not 2 (%zu)\n",
              types.Size());
    return -1;
  }
  Bmap::iterator it= bpmap_.find( types );
  if (it == bpmap_.end()) {
    // New bond parm
    bpmap_.insert( Bpair(types, bp) );
    mprintf("\tAdded new bond params for %s - %s\n", *(types[0]), *(types[1]));
  } else {
    if (allowUpdate) {
      mprintf("\tUpdating bond parameters for %s - %s\n", *(types[0]), *(types[1]));
      it->second = bp;
    } else {
      mprinterr("Error: Update of bond params for %s - %s not allowed.\n",
                *(types[0]), *(types[1]));
      return 1;
    }
  }
  return 0;
}
