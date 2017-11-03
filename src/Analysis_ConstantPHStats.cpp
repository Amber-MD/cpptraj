#include "Analysis_ConstantPHStats.h"
#include "CpptrajStdio.h"
#include "DataSet_PH.h"

// Analysis_ConstantPHStats::Help()
void Analysis_ConstantPHStats::Help() const {

}

// Analysis_ConstantPHStats::Setup()
Analysis::RetType Analysis_ConstantPHStats::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  // Get DataSets
  DataSetList tempDSL;
  std::string dsarg = analyzeArgs.GetStringNext();
  while (!dsarg.empty()) {
    tempDSL += setup.DSL().GetMultipleSets( dsarg );
    dsarg = analyzeArgs.GetStringNext();
  }
  // Remove non-pH data sets
  for (DataSetList::const_iterator ds = tempDSL.begin(); ds != tempDSL.end(); ++ds)
    if ( (*ds)->Type() == DataSet::PH )
      inputSets_.AddCopyOfSet( *ds );
    else
      mprintf("Warning: Set '%s' is not a pH data set, skipping.\n", (*ds)->legend());
  if (inputSets_.empty()) {
    mprinterr("Error: No pH data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    CONSTANT PH STATS:\n");
  inputSets_.List();
  return Analysis::OK;
}

// Analysis_ConstantPHStats::Analyze()
Analysis::RetType Analysis_ConstantPHStats::Analyze() {
  // Loop over all data sets
  for (DataSetList::const_iterator ds = inputSets_.begin(); ds != inputSets_.end(); ++ds)
  {
    DataSet_PH const& PH = static_cast<DataSet_PH const&>( *((DataSet_PH*)*ds) );
    unsigned int nframes = PH.Nframes();
    if (nframes > 0) {
      // Loop over all residues
      for (DataSet_PH::const_iterator res = PH.begin(); res != PH.end(); ++res)
      {
         // Initial state.
         int last_state = res->States().front();
         int n_transitions = 0;
         int n_prot = (int)(res->IsProtonated(last_state));
         int tot_prot = res->Nprotons(last_state);
         // Loop over frames after initial.
         for (unsigned int n = 1; n != res->Nframes(); n++)
         {
           //if ( res->State(n) != last_state )
           if ( res->IsProtonated( res->State(n) ) != res->IsProtonated( last_state ) )
             n_transitions++;
           if ( res->IsProtonated( res->State(n) ) )
             n_prot++;
           tot_prot += res->Nprotons( res->State(n) );
           last_state = res->State(n);
         }
         rprintf("DEBUG: %s '%s %i' n_transitions= %i  n_prot= %i  tot_prot= %i\n",
                 PH.legend(), *(res->Name()), res->Num(), n_transitions, n_prot, tot_prot);
       } // END loop over residues
     }
   } // END loop over DataSets

  return Analysis::OK;
}
