#include "Analysis_ConstantPHStats.h"
#include "CpptrajStdio.h"

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
    if ( (*ds)->Type() == DataSet::PH ) {
      /// Require residue data.
      if ( ((DataSet_PH*)(*ds))->Residues().empty() ) {
        mprinterr("Error: pH set '%s' has no residue info.\n", (*ds)->legend());
        return Analysis::ERR;
      }
      inputSets_.AddCopyOfSet( *ds );
    } else
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
  StatMap Stats;
  // Loop over all data sets
  for (DataSetList::const_iterator ds = inputSets_.begin(); ds != inputSets_.end(); ++ds)
  {
    DataSet_PH const& PH = static_cast<DataSet_PH const&>( *((DataSet_PH*)*ds) );
    if (PH.Nframes() > 0) {
      // Loop over all residues
      for (DataSet_PH::const_iterator res = PH.begin(); res != PH.end(); ++res)
      {
        // Initial state.
        int last_state = res->States().front();
        PHresMap::iterator ph_res;
        // Try to find residue in map.
        StatMap::iterator it = Stats.lower_bound( res->Num() );
        if ( it == Stats.end() || it->first != res->Num() ) {
          // New residue. First create map of solvent pH to residue.
          PHresMap tmp;
          tmp.insert( PHresPair(PH.Solvent_pH(), ResStat(*res, last_state)) );
          it = Stats.insert( it, StatPair(res->Num(), tmp) );
          ph_res = it->second.begin();
        } else {
          // Existing residue. Find pH.
          ph_res = it->second.lower_bound( PH.Solvent_pH() );
          if (ph_res == it->second.end() ||
              ph_res->first != PH.Solvent_pH()) // TODO fix comparison
          {
            // New pH value.
            ph_res = it->second.insert( ph_res, PHresPair(PH.Solvent_pH(),
                                                          ResStat(*res, last_state)) );
          }
        }
        ResStat& stat = ph_res->second;

        // Loop over frames after initial.
        for (unsigned int n = 1; n != res->Nframes(); n++)
        {
          //if ( res->State(n) != last_state )
          if ( res->IsProtonated( res->State(n) ) != res->IsProtonated( last_state ) )
            stat.n_transitions_++;
          if ( res->IsProtonated( res->State(n) ) )
            stat.n_prot_++;
          stat.tot_prot_ += res->Nprotons( res->State(n) );
          last_state = res->State(n);
        }
        rprintf("DEBUG: %s '%s %i' n_transitions= %i  n_prot= %i  tot_prot= %i\n",
                PH.legend(), *(res->Name()), res->Num(),
                stat.n_transitions_, stat.n_prot_, stat.tot_prot_);
      } // END loop over residues
    }
  } // END loop over DataSets

  mprintf("#%-5s %4s %6s %8s %8s %8s\n", "pH", "Name", "Num", "Ntrans", "Nprot", "TotProt");
  int nstats = 0;
  for (StatMap::const_iterator res_map = Stats.begin(); res_map != Stats.end(); ++res_map)
  {
    //int resnum = res_map->first;
    for (PHresMap::const_iterator ph_res = res_map->second.begin();
                                  ph_res != res_map->second.end(); ++ph_res)
    {
      ResStat const& stat = ph_res->second;
      nstats++;
      rprintf("%6.2f %4s %6i %8i %8i %8i\n", ph_res->first,
              *(stat.name_), stat.num_, stat.n_transitions_, stat.n_prot_, stat.tot_prot_);
    }
  }
/*
# ifdef MPI
  // For doing things like pH plots gather all data to a single thread.
  if (!Parallel::EnsembleComm().IsNull() && Parallel::EnsembleComm().Size() > 1)
  {
    if (Parallel::EnsembleComm().Master() {
      std::vector<int> Nstats_on_thread( Parallel::EnsembleComm().Size() );
      Parallel::EnsembleComm().GatherMaster( &nstats, 1, MPI_INT, &Nstats_on_thread[0] );

    } else {
      // Not ensemble master. Send master how many stats I have.
      Parallel::EnsembleComm().GatherMaster( &nstats, 1, MPI_INT, 0 );
      // Send master each stat.
      for (StatMap::const_iterator res_map = Stats.begin(); res_map != Stats.end(); ++res_map)
      {
        //int resnum = res_map->first;
        for (PHresMap::const_iterator ph_res = res_map->second.begin();
                                      ph_res != res_map->second.end(); ++ph_res)
        {
          float ph = 
          ResStat const& stat = ph_res->second;
  */

  // Create a titration curve for each residue.
          

  return Analysis::OK;
}
// =============================================================================
/*
bool Analysis_ConstantPHStats::ResStat::operator==(const ResStat& rhs) const {
  float diff = pH_ - rhs.pH_;
  if (diff < 0.0) diff = -diff;
  return (diff < Constants::SMALL);
}*/
