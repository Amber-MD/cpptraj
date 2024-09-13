#include "EnergyDecomposer.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSetList.h"
#include "../ParameterTypes.h"
#include <algorithm> //std::sort

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

// -----------------------------------------------------------------------------
/** Set up bonds. */
void EnergyDecomposer::setupBonds(BndArrayType const& bondsIn) {
  for (BndArrayType::const_iterator bnd = bondsIn.begin(); bnd != bondsIn.end(); ++bnd)
  {
    if ( selectedAtoms_.AtomInCharMask( bnd->A1() ) ||
         selectedAtoms_.AtomInCharMask( bnd->A2() ) )
    {
      bonds_.push_back( *bnd );
    }
  }
}

/** Topology-based setup.
  * \return 0 if setup OK, 1 if error, -1 if nothing selected.
  */
int EnergyDecomposer::SetupDecomposer(Topology const& topIn) {
  // First set up the mask
  if (topIn.SetupCharMask( selectedAtoms_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", selectedAtoms_.MaskString());
    return 1;
  }
  selectedAtoms_.MaskInfo();
  if (selectedAtoms_.None()) {
    mprintf("Warning: Nothing selected by mask '%s'\n", selectedAtoms_.MaskString());
    return -1;
  }
  // Set up calculation arrays
  if (indices_.empty()) {
    // First time setup
    indices_.reserve( selectedAtoms_.Nselected() );
    for (int idx = 0; idx != topIn.Natom(); idx++)
      if (selectedAtoms_.AtomInCharMask( idx ))
        indices_.push_back( idx );
    energies_.resize( selectedAtoms_.Nselected() );
  } else {
    // Already setup. Warn if indices have changed.
    if ((unsigned int)selectedAtoms_.Nselected() != indices_.size()) {
      // FIXME implement this
      mprinterr("Error: Number of selected atoms has changed in topology '%s'\n", topIn.c_str());
      mprinterr("Error: Not yet supported by energy decomposition.\n");
      return 1;
    }
  }
  // Set up bonds
  bonds_.clear();
  setupBonds( topIn.Bonds() );
  setupBonds( topIn.BondsH() );
  std::sort( bonds_.begin(), bonds_.end() );
  mprintf("DEBUG: Bonds:\n");
  for (BndArrayType::const_iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
    mprintf("\t%s - %s\n", topIn.AtomMaskName(bnd->A1()).c_str(), topIn.AtomMaskName(bnd->A2()).c_str());

  return 0;
}

// -----------------------------------------------------------------------------
