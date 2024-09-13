#include "EnergyDecomposer.h"
#include "../ArgList.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EnergyDecomposer::EnergyDecomposer()
{ }

/** Initialize decomposer. */
int EnergyDecomposer::InitDecomposer(ArgList&, int debugIn)
{
  debug_ = debugIn;

  return 0;
}
