#include <cstdlib> // abs
#include "ForLoop_dataSetBlocks.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "DataSet.h"

/// CONSTRUCTOR
ForLoop_dataSetBlocks::ForLoop_dataSetBlocks() :
  sourceSet_(0),
  blocksize_(0),
  blockoffset_(0),
  idx_(0)
{}

int ForLoop_dataSetBlocks::SetupFor(CpptrajState& State, std::string const& expr, ArgList& argIn)
{
  // <var> datasetblocks <set> blocksize <#> [blockoffset <#>]
  std::string dsname = argIn.GetStringKey("datasetblocks");
  if (dsname.empty()) {
    mprinterr("Error: No data set name given.\n");
    return 1;
  }
  sourceSet_ = State.DSL().GetDataSet( dsname );
  if (sourceSet_ == 0) {
    mprinterr("Error: No data set found with name '%s'\n", dsname.c_str());
    return 1;
  }
  if (sourceSet_->Size() < 1) {
    mprinterr("Error: Data set is empty.\n");
    return 1;
  }
  if (sourceSet_->Group() != DataSet::SCALAR_1D &&
      sourceSet_->Type() != DataSet::VECTOR)
  {
    mprinterr("Error: Set '%s' is not 1D scalar or vector.\n");
    return 1;
  }
  int bs = argIn.getKeyInt("blocksize", 0);
  if (bs < 1) {
    mprinterr("Error: No blocksize or invalid blocksize: %i\n", bs);
    return 1;
  }
  blocksize_ = (unsigned int)bs;
  if (blocksize_ >= sourceSet_->Size()) {
    // TODO make sure we are doing at least 2 iterations?
    mprinterr("Error: block size %u is greater than data set size %zu\n",
              blocksize_, sourceSet_->Size());
    return 1;
  }
  // TODO figure out a clever way to make sure we're not accessing another 'for' block
  blockoffset_ = argIn.getKeyInt("blockoffset", 0);
  if (blockoffset_ == 0) {
    mprintf("Warning: 'blockoffset' not specified, using 'blocksize'.\n");
    blockoffset_ = (int)blocksize_;
  }
  if (abs(blockoffset_) > (int)blocksize_) {
    mprinterr("Error: Block offset is greater than block size.\n");
    return 1;
  }
  // Set up loop variable
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;
  // Description
  std::string description(VarName() + " datasetblocks " + sourceSet_->Meta().Name());
  SetDescription( description );
  return 0;
}
