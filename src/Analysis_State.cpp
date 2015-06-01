#include "Analysis_State.h"
#include "CpptrajStdio.h"

void Analysis_State::Help() {
  mprintf("\t{state <ID>,<dataset>,<min>,<max>} [out <file>]\n"
          "  Data for the specified data set(s) that matches the given criteria\n"
          "  will be assigned state <#>.\n");
}

Analysis::RetType Analysis_State::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  // Get definitions of states if present.
  // Define states as 'state <#>,<dataset>,<min>,<max>'
  std::string state_arg = analyzeArgs.GetStringKey("state");
  if (!state_arg.empty()) {
    while (!state_arg.empty()) {
      // Expect form <#>,<dataset>
      ArgList argtmp(state_arg, ",");
      if (argtmp.Nargs() != 4) {
        mprinterr("Error: Malformed state argument '%s': expect <ID>,<dataset>,<min>,<max>\n",
                  state_arg.c_str());
        return Analysis::ERR;
      }
      std::string state_id = argtmp.GetStringNext();
      // TODO: Check duplicate names
      if (state_id.empty()) {
        mprinterr("Error: In state argument, could not get ID.\n");
        return Analysis::ERR;
      }
      DataSet* ds = datasetlist->GetDataSet( argtmp.GetStringNext() );
      if (ds == 0) return Analysis::ERR;
      if (ds->Ndim() != 1) {
        mprinterr("Error: Only 1D data sets allowed.\n");
        return Analysis::ERR;
      }
      double min = argtmp.getNextDouble(0.0);
      double max = argtmp.getNextDouble(0.0);
      if (max < min) {
        mprinterr("Error: max value cannot be less than min.\n");
        return Analysis::ERR;
      }
      States_.push_back( StateType(state_id, (DataSet_1D*)ds, min, max) );
      state_arg = analyzeArgs.GetStringKey("state");
    }
  }
  if (States_.empty()) {
    mprinterr("Error: No states defined.\n");
    return Analysis::ERR;
  }
  stateOut_ = datasetlist->AddSet(DataSet::INTEGER, analyzeArgs.GetStringNext(), "State");
  if (stateOut_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddSet( stateOut_ );

  mprintf("    STATE: The following states have been set up:\n");
  for (StateArray::const_iterator state = States_.begin(); state != States_.end(); ++state)
    mprintf("\t%u: %20s %12.4f < %-20s < %12.4f\n", state->DS().legend(),
            state - States_.begin(), state->Min(), state->id(), state->Max());
  if (outfile != 0)
    mprintf("\tOutput to file '%s'\n", outfile->DataFilename().full());

  return Analysis::OK;
}

Analysis::RetType Analysis_State::Analyze() {
  // Only process as much data as is in the smallest data set.
  size_t nframes = 0;
  for (StateArray::const_iterator state = States_.begin(); state != States_.end(); ++state)
  {
    if ( state != States_.begin() && nframes != state->DS().Size() )
      mprintf("Warning: Set '%s' for state '%s' has a different # of frames (%zu)"
              " than previous sets (%zu).\n", state->DS().legend(), state->id(),
              state->DS().Size(), nframes);
    nframes = std::min( nframes, state->DS().Size() );
  }
  mprintf("\tProcessing %zu frames.\n", nframes);
  if (nframes < 1) return Analysis::ERR;

  for (size_t frm = 0; frm != nframes; frm++) {
    // Determine which state we are in.
    int state_num = -1;
    for (StateArray::const_iterator state = States_.begin(); state != States_.end(); ++state)
    {
      double dVal = state->DS().Dval( frm );
      if (dVal > state->Min() && dVal < state->Max()) { // TODO: Periodic
        if (state_num != -1)
          mprintf("Warning: Frame %zu already defined as state '%s', also matches state '%s'.\n",
                  frm, States_[state_num].id(), state->id());
        else
          state_num = (int)(state - States_.begin());
      }
    }
    stateOut_->Add( frm, &state_num );
  }

  return Analysis::OK;
} 
