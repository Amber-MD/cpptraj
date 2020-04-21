#include "PotentialTerm_Bond.h"
#include "Topology.h"
#include "CharMask.h"

int PotentialTerm_Bond::SetupTerm(Topology const& topIn, CharMask const& maskIn)
{
  activeBonds_.clear();
  for (BondArray::const_iterator bnd = topIn.Bonds().begin(); bnd != topIn.Bonds().end(); ++bnd)
  {
    if (maskIn.AtomInCharMask( bnd->A1() ) ||
        maskIn.AtomInCharMask( bnd->A2() ))
    {
      //mprintf("DEBUG: Bond %i to %i\n", bnd->A1()+1, bnd->A2()+1);
      activeBonds_.push_back( *bnd );
    }
  }

  return 0;
}
