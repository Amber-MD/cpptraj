#include "ForLoop_dataSetBlocks.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "DataSet.h"
#include "StringRoutines.h"
#ifdef MPI
# include "Parallel.h"
#endif
#include <cstdlib>

/// CONSTRUCTOR
ForLoop_dataSetBlocks::ForLoop_dataSetBlocks() :
  sourceSet_(0),
  currentSet_(0),
  blocksize_(0),
  blockoffset_(0),
  idx_(0),
  mode_(BLOCKS)
{}

void ForLoop_dataSetBlocks::helpText() {
  mprintf("\t<var> datasetblocks <set> blocksize <#> [blockoffset <#>]\n"
          "\t[cumulative [firstblock <#>]]\n"
          "  Loop over blocks of a dataset. The blocks can be either\n"
          "  fixed or cumulative.\n");
}

int ForLoop_dataSetBlocks::SetupFor(CpptrajState& State, ArgList& argIn)
{
  // <var> datasetblocks <set> blocksize <#> [blockoffset <#>] [cumulative [firstblock <#>]]
# ifdef MPI
  if (Parallel::World().Size() > 1) {
    mprinterr("Error: 'datasetblocks' loop cannot be used in parallel.\n");
    return 1;
  }
# endif
  if (argIn.hasKey("cumulative"))
    mode_ = CUMULATIVE;
  else
    mode_ = BLOCKS;
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
  if (mode_ == BLOCKS) {
    if (blockoffset_ == 0) {
      mprintf("Warning: 'blockoffset' not specified, using 'blocksize'.\n");
      blockoffset_ = blocksize_;
    }
  } else if (mode_ == CUMULATIVE) {
    if (blockoffset_ != 0) {
      mprinterr("Error: 'blockoffset' cannot be specified with 'cumulative'\n");
      return 1;
    }
    blockoffset_ = blocksize_;
    blocksize_ = argIn.getKeyInt("firstblock", blocksize_);
  }
  idx_ = argIn.getKeyInt("blockstart", 0);
  // Set up loop variable
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Description
  std::string description("(" + VarName() + " datasetblocks " + sourceSetName_ + ")");
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
      sourceSet_->Group() != DataSet::VECTOR_1D)
  {
    mprinterr("Error: Set '%s' is not 1D scalar or vector.\n", sourceSet_->legend());
    return 1;
  }
  if (sourceSet_->Size() < 1) {
    mprinterr("Error: Set '%s' is empty.\n", sourceSet_->legend());
    return 1;
  }
  // Determine # iterations
  long int total_size = (long int)sourceSet_->Size() - idx_;
  long int niterations = 0;
  if (mode_ == CUMULATIVE) {
    niterations = 1;
    total_size -= blocksize_;
    if (total_size < 0) total_size = 0;
  }
  niterations += ( total_size / labs(blockoffset_) );
  if ( (total_size % labs(blockoffset_)) > 0 )
    ++niterations;
  return niterations;
}

bool ForLoop_dataSetBlocks::EndFor(DataSetList& DSL) {
  // Determine stop of the current block.
  long int block_end = idx_ + blocksize_;
  long int dsidx = 0;
  std::string aspect = "";
  // Check if done.
  switch (mode_) {
    case BLOCKS:
      // Check if done by seeing if the current start value is outside the data set.
      if (idx_ < 0 || idx_ >= (long int)sourceSet_->Size()) return true;
      dsidx = idx_ + 1;
      aspect = "block";
      break;
    case CUMULATIVE:
      // Check if done by seeing if current end value is outside the data set.
      //if (block_end < 0) {
      //  if ( block_end < blockoffset_ ) return true;
      //} else if (block_end > (long int)sourceSet_->size()) {
        //mprintf("DEBUG: Block end %li set size %zu\n", block_end, sourceSet_->Size());
        if ( block_end >= (long int)sourceSet_->Size() + blockoffset_ ) return true;
        // Ensure block end is not too high
        if (block_end > (long int)sourceSet_->Size())
          block_end = (long int)sourceSet_->Size();
        dsidx = block_end;
        aspect = "cumul";
      break;
  }
  //mprintf("DEBUG: Block %li to %li\n", idx_, block_end);
  // Create the subset
  currentSet_ = DSL.AddSet(sourceSet_->Type(), MetaData(VarName(), aspect, dsidx));
  if (currentSet_ == 0) {
    mprinterr("Error: Could not create dataSetBlocks subset.\n");
    return true;
  }
  // Update the set legend
  currentSet_->SetLegend(sourceSet_->Meta().Name() + "_" +
                         integerToString(idx_+1) + "-" +
                         integerToString(block_end));
  // Update the set name in loop variable
  DSL.UpdateStringVar( VarName(), currentSet_->Meta().PrintName() );
  // Allocate memory
  DataSet::SizeArray dsSize(1, blocksize_);
  if ( currentSet_->MemAlloc( dsSize ) ) {
    mprinterr("Internal Error: Could not allocate dataSetBlocks subset; not yet supported.\n");
    return true;
  }
  // Copy block
  currentSet_->CopyBlock(0, sourceSet_, idx_, blocksize_);
  switch (mode_) {
    case BLOCKS:
      // Increment the start
      idx_ += blockoffset_;
      break;
    case CUMULATIVE:
      // Increment the size
      blocksize_ += blockoffset_;
      break;
  }
  return false;
}
