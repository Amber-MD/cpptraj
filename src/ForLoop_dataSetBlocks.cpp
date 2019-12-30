#include "ForLoop_dataSetBlocks.h"
#include "CpptrajStdio.h"

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
  blocksize_ = argIn.getKeyInteger("blocksize", 0);
  if (blocksize_ < 1) {
    mprinterr("Error: No blocksize or invalid blocksize: %i\n", blocksize_);
    return 1;
  }
  if ((unsigned int)blocksize_ >= sourceSet_->Size()) {
    // TODO make sure we are doing at least 2 iterations?
    mprinterr("Error: block size %i is greater than data set size %zu\n",
              blocksize_, sourceSet_=>Size());
    return 1;
  }
  // TODO figure out a clever way to make sure we're not accessing another 'for' block
  blockoffset_ = argIn.getKeyInteger("blockoffset", -1);
  if (blockoffset < 1)
    block
