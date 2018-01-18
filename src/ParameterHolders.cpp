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
  // Check if bond parm for these types exist
  Bmap::iterator it = bpmap_.begin();
  for (; it != bpmap_.end(); ++it)
    if (it->first == types) break;
  if (it == bpmap_.end()) {
    // New bond parm
    bpmap_.push_back( Bpair(types, bp) );
    mprintf("\tAdded new bond params for %s - %s\n", *(types[0]), *(types[1]));
  } else {
    if (allowUpdate) {
      mprintf("\tUpdating bond parameters for %s - %s\n", *(types[0]), *(types[1]));
      it->second = bp;
    } else {
      mprinterr("Error: Update of bond params for %s - %s not allowed.\n",
                *(types[0]), *(types[1]));
      mprinterr("DBG: %s %s\n", *(it->first[0]), *(it->first[1]));
      if (types == it->first) mprinterr("DBG: MATCH.\n");
      //if (types < it->first) mprinterr("DBG: LESS THAN.\n");
      return 1;
    }
  }
  return 0;
}
