#include <cmath>
#include <algorithm>
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
  // Loop over all data sets
  for (DataSetList::const_iterator ds = inputSets_.begin(); ds != inputSets_.end(); ++ds)
  {
    DataSet_pH* phset = (DataSet_pH*)*ds;
    if (phset->Size() > 0) {
      Stats_.push_back( ResStat(phset, phset->State(0)) );
      DataSet_pH const& PH = static_cast<DataSet_pH const&>( *phset );
      // Get residue stats
      ResStat& stat = Stats_.back(); 
      // Initial state.
      int last_state = PH.State(0);
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

  // Print some results grouped by residue
  std::sort( Stats_.begin(), Stats_.end(), num_ph_sort() );
  mprintf("#%-5s %4s %6s %8s %8s %8s\n", "pH", "Name", "Num", "Ntrans", "Nprot", "TotProt");
  for (Rarray::const_iterator stat = Stats_.begin(); stat != Stats_.end(); ++stat)
  {
    DataSet_pH const& PH = static_cast<DataSet_pH const&>( *(stat->ds_) );
    rprintf("%6.2f %4s %6i %8i %8i %8i\n", PH.Solvent_pH(), *(PH.Res().Name()), PH.Res().Num(),
            stat->n_transitions_, stat->n_prot_, stat->tot_prot_);
  }

  // Print cphstats-like output, grouped by pH
  if (statsOut_ != 0) {
    std::sort( Stats_.begin(), Stats_.end(), ph_num_sort() );
    float current_pH = -100.0;
    bool write_pH = true;
    unsigned int tot_prot = 0;
    for (Rarray::const_iterator stat = Stats_.begin(); stat != Stats_.end(); ++stat)
    {
      if (write_pH) {
        current_pH = stat->ds_->Solvent_pH();
        statsOut_->Printf("Solvent pH is %8.3f\n", current_pH);
        write_pH = false;
        tot_prot = 0;
      }
      statsOut_->Printf("%3s %-4i", stat->ds_->Res().Name().Truncated().c_str(),
                        stat->ds_->Res().Num());
      double dsize = (double)stat->ds_->Size();
      double dnprot = (double)stat->n_prot_;
      if (dnprot > 0.0) {
        double pKa = (double)current_pH - log10( (dsize - dnprot) / dnprot );
        double offset = pKa - current_pH;
        statsOut_->Printf(": Offset %6.3f", offset);
        statsOut_->Printf("  Pred %6.3f", pKa);
      } else {
        statsOut_->Printf(": Offset %6s", "inf");
        statsOut_->Printf("  Pred %6s", "inf");
      }
      statsOut_->Printf("  Frac Prot %5.3f", dnprot / dsize);
      statsOut_->Printf(" Transitions %9i\n", stat->n_transitions_);
      tot_prot += stat->tot_prot_;

      Rarray::const_iterator next = stat + 1;
      if (next == Stats_.end() || next->ds_->Solvent_pH() > current_pH) {
        write_pH = true;
        statsOut_->Printf("\nAverage total molecular protonation: %7.3f\n",
                         (double)tot_prot / dsize);
      }
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
