#include "EnergyDecomposer.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSetList.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EnergyDecomposer::EnergyDecomposer() :
  eneOut_(0)
{ }

/** Initialize decomposer. */
int EnergyDecomposer::InitDecomposer(ArgList& argIn, DataSetList& DSLin, int debugIn)
{
  debug_ = debugIn;

  // Get atom mask
  if (selectedAtoms_.SetMaskString( argIn.GetMaskNext() ))
    return 1;
  // Output DataSet
  std::string setname = argIn.GetStringNext();
  if (setname.empty())
    setname = DSLin.GenerateDefaultName("DECOMP");
  eneOut_ = DSLin.AddSet(DataSet::XYMESH, MetaData(setname));
  if (eneOut_ == 0) {
    mprinterr("Error: Could not allocate decomp. output set '%s'\n", setname.c_str());
    return 1;
  }

  return 0;
}

/** Print options to stdout. */
void EnergyDecomposer::PrintOpts() const {
  if (eneOut_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::PrintOpts() called before initialization.\n");
    return;
  }
  mprintf("\tCalculating for atoms selected by mask: %s\n", selectedAtoms_.MaskString());
  mprintf("\tData set name: %s\n", eneOut_->legend());
}
