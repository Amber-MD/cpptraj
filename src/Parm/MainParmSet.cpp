#include "MainParmSet.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Parameters.h"
#include "../DataSetList.h"

using namespace Cpptraj::Parm;

/** CONSTRUCTOR */
MainParmSet::MainParmSet() :
  mainParmSet_(0),
  free_parmset_mem_(false)
{}

/** DESTRUCTOR */
MainParmSet::~MainParmSet() {
  if (mainParmSet_ != 0 && free_parmset_mem_)
    delete mainParmSet_;
}

const char* MainParmSet::parm_keywords_ = "parmset <parameter setname>";

/** Initialize */
int MainParmSet::InitMainParmSet(ArgList& argIn, DataSetList const& DSL, int debugIn)
{
  t_total_.Start();
  debug_ = debugIn;

  if (getParameterSets(argIn, DSL)) return 1;

  t_total_.Stop();
  return 0;
}

/** Write timing info to stdout. */
void MainParmSet::TimingInfo(double total, int indent) const {
  t_total_.WriteTiming(indent, "Get Parameters:", total);
}

/** Get parameter sets. */
int MainParmSet::getParameterSets(ArgList& argIn, DataSetList const& DSL) {
  // Clear any existing set
  if (mainParmSet_ != 0) {
    if (free_parmset_mem_) delete mainParmSet_;
  }
  mainParmSet_ = 0;
  free_parmset_mem_ = false;
  // Look for parmset args
  typedef std::vector<DataSet_Parameters*> Parray;
  Parray ParamSets;
  std::string parmset = argIn.GetStringKey("parmset");
  if (parmset.empty()) {
    mprintf("\tNo parameter set(s) specified with 'parmset'; using any loaded sets.\n");
    // See if there are any parameter sets.
    DataSetList sets = DSL.GetSetsOfType( "*", DataSet::PARAMETERS );
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
      ParamSets.push_back( (DataSet_Parameters*)(*it) );
  } else {
    while (!parmset.empty()) {
      DataSetList sets = DSL.GetSetsOfType( parmset, DataSet::PARAMETERS );
      if (sets.empty()) {
        mprintf("Warning: No parameter sets corresponding to '%s'\n", parmset.c_str());
      } else {
        for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
          ParamSets.push_back( (DataSet_Parameters*)(*it) );
      }
      parmset = argIn.GetStringKey("parmset");
    }
  }
  //if (ParamSets.empty()) {
  //  mprinterr("Error: No parameter sets.\n");
  //  return CpptrajState::ERR;
  //}
  if (!ParamSets.empty()) {
    mprintf("\tParameter sets:\n"); // TODO put these names in the final combined parameter set
    for (Parray::const_iterator it = ParamSets.begin(); it != ParamSets.end(); ++it)
      mprintf("\t  %s\n", (*it)->legend());

    // Combine parameters if needed

    if (ParamSets.size() == 1)
      mainParmSet_ = ParamSets.front();
    else {
      free_parmset_mem_ = true;
      mprintf("\tCombining parameter sets.\n");
      Parray::const_iterator it = ParamSets.begin();
      mainParmSet_ = new DataSet_Parameters( *(*it) );
      ++it;
      Cpptraj::Parm::ParameterSet::UpdateCount UC;
      for (; it != ParamSets.end(); ++it)
        mainParmSet_->UpdateParamSet( *(*it), UC, debug_, debug_+1 ); // Make it so verbosity is at least 1 to report overwritten params 
    }
  }
  return 0;
}
