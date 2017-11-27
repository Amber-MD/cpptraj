#include "Analysis_ConstantPHStats.h"
#include "CpptrajStdio.h"
#include "DataSet_pH.h"

// Analysis_ConstantPHStats::Help()
void Analysis_ConstantPHStats::Help() const {
  mprintf("\t<pH sets> [statsout <statsfile>]\n");
}

// Analysis_ConstantPHStats::Setup()
Analysis::RetType Analysis_ConstantPHStats::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  statsOut_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("statsout"), 
                                         "Constant pH stats", DataFileList::TEXT);
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
      if ( ((DataSet_pH*)(*ds))->Res().Num() == -1 ) {
        mprinterr("Error: pH set '%s' has no residue info.\n", (*ds)->legend());
        return Analysis::ERR;
      }
      inputSets_.AddCopyOfSet( *ds );
    } else if ( (*ds)->Type() == DataSet::PH_REMD ) {
      mprinterr("Error: pH set '%s' must be sorted first.\n", (*ds)->legend());
      return Analysis::ERR;
    } else
      mprintf("Warning: Set '%s' is not a pH data set, skipping.\n", (*ds)->legend());
  if (inputSets_.empty()) {
    mprinterr("Error: No pH data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    CONSTANT PH STATS:\n");
  if (statsOut_ != 0)
    mprintf("\tConstant pH statistics (cphstats style) output to '%s'\n",
            statsOut_->Filename().full());
  mprintf("\tInput pH data sets:\n");
  inputSets_.List();
  return Analysis::OK;
}

// Analysis_ConstantPHStats::Analyze()
Analysis::RetType Analysis_ConstantPHStats::Analyze() {
  StatMap Stats;
  // Loop over all data sets
  for (DataSetList::const_iterator ds = inputSets_.begin(); ds != inputSets_.end(); ++ds)
  {
    DataSet_pH const& PH = static_cast<DataSet_pH const&>( *((DataSet_pH*)*ds) );
    if (PH.Size() > 0) {
        // Initial state.
        int last_state = PH.State(0);
        PHresMap::iterator ph_res;
        // Try to find residue in map.
        StatMap::iterator it = Stats.lower_bound( PH.Res().Num() );
        if ( it == Stats.end() || it->first != PH.Res().Num() ) {
          // New residue. First create map of solvent pH to residue.
          PHresMap tmp;
          tmp.insert( PHresPair(PH.Solvent_pH(), ResStat(PH.Res(), last_state)) );
          it = Stats.insert( it, StatPair(PH.Res().Num(), tmp) );
          ph_res = it->second.begin();
        } else {
          // Existing residue. Find pH.
          ph_res = it->second.lower_bound( PH.Solvent_pH() );
          if (ph_res == it->second.end() ||
              ph_res->first != PH.Solvent_pH()) // TODO fix comparison
          {
            // New pH value.
            ph_res = it->second.insert( ph_res, PHresPair(PH.Solvent_pH(),
                                                          ResStat(PH.Res(), last_state)) );
          }
        }
        ResStat& stat = ph_res->second;

        // Loop over frames after initial.
        for (unsigned int n = 1; n != PH.Size(); n++)
        {
          //if ( PH.State(n) != last_state )
          if ( PH.Res().IsProtonated( PH.State(n) ) != PH.Res().IsProtonated( last_state ) )
            stat.n_transitions_++;
          if ( PH.Res().IsProtonated( PH.State(n) ) )
            stat.n_prot_++;
          stat.tot_prot_ += PH.Res().Nprotons( PH.State(n) );
          last_state = PH.State(n);
        }
        rprintf("DEBUG: %s '%s %i' n_transitions= %i  n_prot= %i  tot_prot= %i\n",
                PH.legend(), *(PH.Res().Name()), PH.Res().Num(),
                stat.n_transitions_, stat.n_prot_, stat.tot_prot_);
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
