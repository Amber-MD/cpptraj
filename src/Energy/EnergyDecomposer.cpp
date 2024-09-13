#include "EnergyDecomposer.h"
#include <cmath> // sqrt for Ene_Bond etc
#include "Ene_Bond.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataFileList.h"
#include "../DataSet_Mesh.h"
#include "../DataSetList.h"
#include "../ParameterTypes.h"
#include <algorithm> //std::sort

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
EnergyDecomposer::EnergyDecomposer() :
  eneOut_(0),
  outfile_(0),
  debug_(0),
  currentTop_(0)
{ }

/** Initialize decomposer. */
int EnergyDecomposer::InitDecomposer(ArgList& argIn, DataSetList& DSLin, DataFileList& DFLin,
                                     int debugIn)
{
  debug_ = debugIn;
  // Process keywords
  outfile_ = DFLin.AddDataFile( argIn.GetStringKey("out"), argIn );
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
# ifdef MPI
  eneOut_->SetNeedsSync( false ); // Not a time series
# endif
  eneOut_->ModifyDim(Dimension::X).SetLabel("Atom");
  if (outfile_ != 0)
    outfile_->AddDataSet( eneOut_ );

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
  if (outfile_ != 0)
    mprintf("\tOutput file: %s\n", outfile_->DataFilename().full());
}

// -----------------------------------------------------------------------------
/** Set up bonds. */
int EnergyDecomposer::setupBonds(BndArrayType const& bondsIn) {
  for (BndArrayType::const_iterator bnd = bondsIn.begin(); bnd != bondsIn.end(); ++bnd)
  {
    if ( selectedAtoms_.AtomInCharMask( bnd->A1() ) ||
         selectedAtoms_.AtomInCharMask( bnd->A2() ) )
    {
      if (bnd->Idx() < 0) {
        mprinterr("Error: Bond %i - %i does not have parameters, cannot calculate energy.\n",
                  bnd->A1()+1, bnd->A2()+1);
        return 1;
      }
      bonds_.push_back( *bnd );
    }
  }
  return 0;
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
  if (energies_.empty()) {
    // First time setup
//    indices_.reserve( selectedAtoms_.Nselected() );
//    for (int idx = 0; idx != topIn.Natom(); idx++)
//      if (selectedAtoms_.AtomInCharMask( idx ))
//        indices_.push_back( idx );
    energies_.resize( topIn.Natom() );
  } else {
    // Already setup. Warn if indices have changed.
    //if ((unsigned int)selectedAtoms_.Nselected() != indices_.size())
    if ((unsigned int)topIn.Natom() != energies_.size())
    {
      // FIXME implement this
      mprinterr("Error: Number of atoms has changed in topology '%s'\n", topIn.c_str());
      mprinterr("Error: Now %i atoms, expected %zu\n", topIn.Natom(), energies_.size());
      mprinterr("Error: Not yet supported by energy decomposition.\n");
      return 1;
    }
  }
  // Set up bonds
  bonds_.clear();
  if (setupBonds( topIn.Bonds() )) return 1;
  if (setupBonds( topIn.BondsH() )) return 1;
  std::sort( bonds_.begin(), bonds_.end() );

  // DEBUG
  mprintf("DEBUG: Saving energy for atoms:\n");
  for (int idx = 0; idx != topIn.Natom(); idx++)
    if (selectedAtoms_.AtomInCharMask( idx ))
      mprintf("\t%s\n", topIn.AtomMaskName( idx ).c_str());
  mprintf("DEBUG: Bonds:\n");
  for (BndArrayType::const_iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
    mprintf("\t%s - %s\n", topIn.AtomMaskName(bnd->A1()).c_str(), topIn.AtomMaskName(bnd->A2()).c_str());

  currentTop_ = &topIn;

  return 0;
}

// -----------------------------------------------------------------------------
/** Save energy contribution to atom if it is selected. */
void EnergyDecomposer::saveEne(int idx, double ene_cont) {
  if (selectedAtoms_.AtomInCharMask( idx ))
    energies_[ idx ].accumulate( ene_cont );
}

/** Calculate bond energies. */
void EnergyDecomposer::calcBonds( Frame const& frameIn ) {
  for (BndArrayType::const_iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
  {
    BondParmType const& BP = currentTop_->BondParm()[ bnd->Idx() ];
    double ene = Ene_Bond<double>( frameIn.XYZ( bnd->A1() ),
                                   frameIn.XYZ( bnd->A2() ),
                                   BP.Req(), BP.Rk() );
    // Divide the energy equally between the two atoms.
    double ene_half = ene * 0.5;
    saveEne( bnd->A1(), ene_half );
    saveEne( bnd->A2(), ene_half );
  }
}

/** Calculate and decompose energies. */
int EnergyDecomposer::CalcEne(Frame const& frameIn) {
  if (currentTop_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::CalcEne() called before setup.\n");
    return 1;
  }
  // Bonds
  calcBonds(frameIn);

  return 0;
}

// -----------------------------------------------------------------------------
/** Finish the calculation by putting the results into the output DataSet. */
int EnergyDecomposer::FinishCalc() {
  if (energies_.empty() || eneOut_ == 0) {
    mprinterr("Internal Error: EnergyDecomposer::FinishCalc() called before setup.\n");
    return 1;
  }
  // Only add entities that have data.
  DataSet_Mesh& set = static_cast<DataSet_Mesh&>( *eneOut_ );
  set.Clear();
  // TODO allocate?
  for (unsigned int idx = 0; idx != energies_.size(); idx++) {
    if ( energies_[idx].nData() > 0 ) {
      set.AddXY( idx+1, energies_[idx].mean() );
    }
  }
  return 0;
}
