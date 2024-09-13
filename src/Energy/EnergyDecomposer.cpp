#include "EnergyDecomposer.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EnergyDecomposer::EnergyDecomposer()
{ }

/** Initialize decomposer. */
int EnergyDecomposer::InitDecomposer(ArgList& argIn, int debugIn)
{
  debug_ = debugIn;

  // Get atom mask
  if (selectedAtoms_.SetMaskString( argIn.GetMaskNext() ))
    return 1;

  return 0;
}

/** Print options to stdout. */
void EnergyDecomposer::PrintOpts() const {
  mprintf("\tCalculating for atoms selected by mask: %s\n", selectedAtoms_.MaskString());
}
