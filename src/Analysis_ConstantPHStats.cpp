#include <cmath>
#include <algorithm>
#include "Analysis_ConstantPHStats.h"
#include "CpptrajStdio.h"
#include "DataSet_pH.h"
#include "DataSet_Mesh.h"
#include "StringRoutines.h"

// Analysis_ConstantPHStats::Help()
void Analysis_ConstantPHStats::Help() const {
  mprintf("\t<pH sets> [statsout <statsfile>] [deprot]\n"
          "\t[fracplot [fracplotout <file>]]\n");
}

// Analysis_ConstantPHStats::Setup()
Analysis::RetType Analysis_ConstantPHStats::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  statsOut_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("statsout"), 
                                         "Constant pH stats", DataFileList::TEXT, false
#                                        ifdef MPI
                                         , Parallel::MasterComm()
#                                        endif
                                        );
  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("CPH");
  createFracPlot_ = analyzeArgs.hasKey("fracplot");
  useFracProtonated_ = !analyzeArgs.hasKey("deprot");
  const char* FRACSTR_;
  if (useFracProtonated_)
    FRACSTR_ = "protonated";
  else
    FRACSTR_ = "deprotonated";
  if (createFracPlot_) {
    fracPlotOut_ = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("fracplotout") );
    if (fracPlotOut_ != 0) {
      fracPlotOut_->ProcessArgs("xlabel pH ylabel \"Frac. " + std::string(FRACSTR_) + "\" noensextension");
#     ifdef MPI
      // Fraction plot should only ever be written by the overall master
      // since it needs data from every ensemble member.
      fracPlotOut_->SetThreadCanWrite( Parallel::MasterComm().Master() );
#     endif
    }
  }
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
    } else if ( (*ds)->Type() == DataSet::PH_EXPL ) {
      mprinterr("Error: pH set '%s' must be sorted first.\n", (*ds)->legend());
      return Analysis::ERR;
    } else
      mprintf("Warning: Set '%s' is not a pH data set, skipping.\n", (*ds)->legend());
  if (inputSets_.empty()) {
    mprinterr("Error: No pH data sets.\n");
    return Analysis::ERR;
  }
  masterDSL_ = setup.DslPtr();

  mprintf("    CONSTANT PH STATS:\n");
  if (statsOut_ != 0)
    mprintf("\tConstant pH statistics (cphstats style) output to '%s'\n",
            statsOut_->Filename().full());
  if (createFracPlot_) {
    mprintf("\tFraction %s vs pH will be calculated", FRACSTR_);
    if (fracPlotOut_ != 0)
      mprintf(" and written to '%s'", fracPlotOut_->DataFilename().full());
    mprintf("\n");
  }
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
    if (phset->Size() > 1) {
      //Stats_.push_back( ResStat(phset, phset->State(0)) );
      // NOTE: Ignore the first state since no dynamics have been done yet.
      Stats_.push_back( ResStat(phset) );
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
      stat.nframes_ += (PH.Size() - 1);
      if (debug_ > 0)
        rprintf("DEBUG: %s '%s %i' n_transitions= %u  n_prot= %u  tot_prot= %u  nframes= %u\n",
                PH.legend(), *(PH.Res().Name()), PH.Res().Num(),
                stat.n_transitions_, stat.n_prot_, stat.tot_prot_, stat.nframes_);
    }
  } // END loop over DataSets

  // Print some results grouped by residue
  std::sort( Stats_.begin(), Stats_.end(), num_ph_sort() );
# ifdef MPI
  DataSetList setsToSync;
  bool needToSync = false;
# endif
  int lastRes = -1;
  DataSet_Mesh* fracPlot = 0;
  mprintf("#%-5s %4s %6s %8s %8s %8s\n", "pH", "Name", "Num", "Ntrans", "Nprot", "TotProt");
  for (Rarray::const_iterator stat = Stats_.begin(); stat != Stats_.end(); ++stat)
  {
    DataSet_pH const& PH = static_cast<DataSet_pH const&>( *(stat->ds_) );
    if (lastRes != PH.Res().Num()) { // Assuming residue numbers are unique
      lastRes = PH.Res().Num();
      // New residue
      if (createFracPlot_) {
        fracPlot = (DataSet_Mesh*)masterDSL_->AddSet(DataSet::XYMESH,
                                                     MetaData(dsname_, "Frac", PH.Res().Num()));
        if (fracPlot == 0) return Analysis::ERR;
        fracPlot->SetLegend(PH.Res().Name().Truncated() + ":" + integerToString(PH.Res().Num()));
#       ifdef MPI
        if (fracPlotOut_ != 0) {
          if (!fracPlotOut_->EnsExt()) {
            needToSync = true;
            fracPlotOut_->AddDataSet(fracPlot);
          } else
            fracPlotOut_->AddDataSet(fracPlot);
        }
        if (needToSync) {
          fracPlot->SetNeedsSync( true );
          setsToSync.AddCopyOfSet( fracPlot );
        }
#       else
        if (fracPlotOut_ != 0) fracPlotOut_->AddDataSet(fracPlot);
#       endif
      }
    }

    rprintf("%6.2f %4s %6i %8u %8u %8u\n", PH.Solvent_pH(), *(PH.Res().Name()), PH.Res().Num(),
            stat->n_transitions_, stat->n_prot_, stat->tot_prot_);

    if (fracPlot != 0) {
      double frac = (double)stat->n_prot_ / (double)stat->nframes_;
      if (!useFracProtonated_)
        frac = 1.0 - frac;
      fracPlot->AddXY( PH.Solvent_pH(), frac );
    }
  }
# ifdef MPI
  if (createFracPlot_)
    setsToSync.SynchronizeData( Parallel::MasterComm() );
# endif

  // Print cphstats-like output, grouped by pH
  if (statsOut_ != 0) {
    std::sort( Stats_.begin(), Stats_.end(), ph_num_sort() );
    float current_pH = -100.0;
    bool write_pH = true;
    unsigned int tot_prot = 0;
#   ifdef MPI
    // If shared, each rank dumps its stats in turn. Otherwise writing
    // to separate files.
    int startRank, stopRank;
    if (statsOut_->IsMPI()) {
      startRank = 0;
      stopRank = Parallel::MasterComm().Size();
    } else {
      startRank = Parallel::MasterComm().Rank();
      stopRank = startRank + 1;
    }
    for (int rank = startRank; rank != stopRank; ++rank) {
      if (rank == Parallel::MasterComm().Rank()) {
#   endif
        for (Rarray::const_iterator stat = Stats_.begin(); stat != Stats_.end(); ++stat)
        {
          if (write_pH) {
            current_pH = stat->ds_->Solvent_pH();
            statsOut_->Printf("Solvent pH is %8.3f\n", current_pH);
            write_pH = false;
            tot_prot = 0;
          }
          statsOut_->Printf("%3s %-4i", *(stat->ds_->Res().Name()),
                            stat->ds_->Res().Num());
          double dsize = (double)stat->nframes_;
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
        } // END loop over Stats_
#   ifdef MPI
      }
      // Hold up before next rank writes. 
      Parallel::MasterComm().Barrier();
    } // END loop over ranks
#   endif
  } // END if statsOut_ != 0
    
  return Analysis::OK;
}
