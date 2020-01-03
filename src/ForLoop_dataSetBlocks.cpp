#include "ForLoop_dataSetBlocks.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "DataSet.h"

/// CONSTRUCTOR
ForLoop_dataSetBlocks::ForLoop_dataSetBlocks() :
  sourceSet_(0),
  currentSet_(0),
  blocksize_(0),
  blockoffset_(0),
  idx_(0)
{}

int ForLoop_dataSetBlocks::SetupFor(CpptrajState& State, ArgList& argIn)
{
  // <var> datasetblocks <set> blocksize <#> [blockoffset <#>]
  sourceSetName_ = argIn.GetStringKey("datasetblocks");
  if (sourceSetName_.empty()) {
    mprinterr("Error: No data set name given.\n");
    return 1;
  }
  blocksize_ = argIn.getKeyInt("blocksize", 0);
  if (blocksize_ < 1) {
    mprinterr("Error: No blocksize or invalid blocksize: %li\n", blocksize_);
    return 1;
  }
  blockoffset_ = argIn.getKeyInt("blockoffset", 0);
  if (blockoffset_ == 0) {
    mprintf("Warning: 'blockoffset' not specified, using 'blocksize'.\n");
    blockoffset_ = blocksize_;
  }
  idx_ = argIn.getKeyInt("blockstart", 0);
  // Set up loop variable
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Description
  std::string description(VarName() + " datasetblocks " + sourceSetName_);
  SetDescription( description );
  return 0;
}

int ForLoop_dataSetBlocks::BeginFor(DataSetList const& DSL) {
  sourceSet_ = DSL.GetDataSet( sourceSetName_ );
  if (sourceSet_ == 0) {
    mprinterr("Error: No data set found with name '%s'\n", sourceSetName_.c_str());
    return 1;
  }
  if (sourceSet_->Group() != DataSet::SCALAR_1D &&
      sourceSet_->Type() != DataSet::VECTOR)
  {
    mprinterr("Error: Set '%s' is not 1D scalar or vector.\n", sourceSet_->legend());
    return 1;
  }
  if (sourceSet_->Size() < 1) {
    mprinterr("Error: Set '%s' is empty.\n", sourceSet_->legend());
    return 1;
  }
  // Determine # iterations
  long int niterations = (long int)sourceSet_->Size() / labs(blockoffset_);
  if ( ((long int)sourceSet_->Size() % labs(blockoffset_)) > 0 )
    ++niterations;
  return niterations;
}

bool ForLoop_dataSetBlocks::EndFor(DataSetList& DSL) {
  // Check if done by seeing if the current start value is outside the data set.
  if (idx_ < 0 || idx_ >= (long int)sourceSet_->Size()) return true;
  // Determine stop of the current block.
  long int block_end = idx_ + blocksize_;
  mprintf("DEBUG: Block %li to %li\n", idx_, block_end);
  // Create the subset
  currentSet_ = DSL.AddSet(sourceSet_->Type(), MetaData(VarName(), idx_));
  if (currentSet_ == 0) {
    mprinterr("Error: Could not create dataSetBlocks subset.\n");
    return true;
  }
  // Update the set name
  DSL.UpdateStringVar( VarName(), currentSet_->Meta().PrintName() );
  // Allocate memory
  DataSet::SizeArray dsSize(1, blocksize_);
  if ( currentSet_->MemAlloc( dsSize ) ) {
    mprinterr("Error: Could not allocate dataSetBlocks subset.\n");
    return true;
  }
  // Copy block
  currentSet_->CopyBlock(0, sourceSet_, idx_, blocksize_);
  // Increment the start
  idx_ += blockoffset_;
  return false;
}
