#include "Action_FilterByData.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

void Action_FilterByData::Help() const {
  mprintf("\t{<dataset arg> min <min> max <max> ...} [out <file>] [name <setname>]\n"
          "\t[multi]\n"
          "  This action has two modes. In the first mode, for all following actions\n"
          "  only frames that are between <min> and <max> of all data sets selected\n"
          "  by each <dataset arg> are allowed to pass. A data set with name <setname>\n"
          "  will be created containing a 1 if the frame passed and 0 if the frame was\n"
          "  filtered out.\n"
          "  If 'multi' is specified then only filter data sets will be created for each\n"
          "  data set instead.\n"
          "  There must be at least one <min> and <max> argument, and can be as many as\n"
          "  there are specified data sets.\n");
}

// Action_FilterByData::Init()
Action::RetType Action_FilterByData::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  Npassed_ = 0;
  Nfiltered_ = 0;
  multi_ = actionArgs.hasKey("multi");
  std::string dsname = actionArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = init.DSL().GenerateDefaultName("Filter");
  DataFile* maxminfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Get min and max args.
  while (actionArgs.Contains("min"))
    Min_.push_back( actionArgs.getKeyDouble("min", 0.0) );
  while (actionArgs.Contains("max"))
    Max_.push_back( actionArgs.getKeyDouble("max", 0.0) );
  if (Min_.empty()) {
    mprinterr("Error: At least one 'min' arg must be specified.\n");
    return Action::ERR;
  }
  if (Max_.empty()) {
    mprinterr("Error: At least one 'max' arg must be specified.\n");
    return Action::ERR;
  }
  if (Min_.size() != Max_.size()) {
    mprinterr("Error: # of 'min' args (%zu) != # of 'max' args (%zu)\n",
              Min_.size(), Max_.size());
    return Action::ERR;
  }
  // Get DataSets from remaining arguments
  Dsets_.AddSetsFromArgs( actionArgs.RemainingArgs(), init.DSL() );
  if (Dsets_.empty()) {
    mprinterr("Error: No data sets specified.\n");
    return Action::ERR;
  }
  // Check number of DataSets and min/max args.
  if ( Dsets_.size() < Min_.size() ) {
    mprinterr("Error: More 'min'/'max' args (%zu) than data sets (%zu).\n",
              Min_.size(), Dsets_.size());
    return Action::ERR;
  }
  if ( Dsets_.size() > Min_.size() ) {
    unsigned int Nremaining = Dsets_.size() - Min_.size();
    double useMin = Min_.back();
    double useMax = Max_.back();
    mprintf("Warning: More data sets than 'min'/'max' args.\n"
            "Warning:  Using min=%f and max=%f for last %u data sets.\n",
            useMin, useMax, Nremaining);
    for (unsigned int ds = 0; ds < Nremaining; ++ds) {
      Min_.push_back( useMin );
      Max_.push_back( useMax );
    }
  }
  // Set up output data set(s)
  if (!multi_) {
    maxmin_ = init.DSL().AddSet( DataSet::INTEGER, dsname );
    if (maxmin_ == 0) return Action::ERR;
    if (maxminfile != 0)
      maxminfile->AddDataSet( maxmin_ );
  } else {
    for (unsigned int idx = 0; idx < Dsets_.size(); idx++) {
      DataSet* ds = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsname, idx));
      if (ds == 0) return Action::ERR;
      ds->SetLegend("Filter("+ds->Meta().PrintName()+")");
      outsets_.push_back( ds );
      if (maxminfile != 0)
        maxminfile->AddDataSet( ds );
    }
  }

  mprintf("    FILTER:");
  if (!multi_)
    mprintf(" Filtering out frames using %zu data sets.\n", Dsets_.size());
  else
    mprintf(" Creating filter data sets for %zu data sets.\n", Dsets_.size());
  for (unsigned int ds = 0; ds < Dsets_.size(); ds++)
    mprintf("\t%.4f < '%s' < %.4f\n", Min_[ds], Dsets_[ds]->legend(), Max_[ds]);
  if (maxminfile != 0)
    mprintf("\tFilter frame info will be written to %s\n", maxminfile->DataFilename().full());
# ifdef MPI
  if (!multi_ && init.TrajComm().Size() > 1)
    mprintf("Warning: Trajectories written after 'filter' may have issues if\n"
            "Warning:   the number of processes writing is > 1 (currently %i processes)\n",
            init.TrajComm().Size());
# endif
  return Action::OK;
}

// Action_FilterByData::DoAction()
Action::RetType Action_FilterByData::DoAction(int frameNum, ActionFrame& frm)
{
  static int ONE = 1;
  static int ZERO = 0;
  if (!multi_) {
    // Check if frame is within max/min of every data set.
    for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    {
      double dVal = Dsets_[ds]->Dval(frm.TrajoutNum());
      //mprintf("DBG: maxmin[%u]: dVal = %f, min = %f, max = %f\n",ds,dVal,Min_[ds],Max_[ds]);
      // If value from dataset not within min/max, exit now.
      if (dVal < Min_[ds] || dVal > Max_[ds]) {
        maxmin_->Add( frameNum, &ZERO );
        Nfiltered_++;
        return Action::SUPPRESS_COORD_OUTPUT;
      }
    }
    maxmin_->Add( frameNum, &ONE );
    Npassed_++;
  } else {
    // For each data set 1 if within min/max, 0 otherwise.
    for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    {
      double dVal = Dsets_[ds]->Dval(frm.TrajoutNum());
      if (dVal < Min_[ds] || dVal > Max_[ds]) 
        outsets_[ds]->Add( frameNum, &ZERO );
      else
        outsets_[ds]->Add( frameNum, &ONE  );
    }
  }
  return Action::OK;
}

/** \return Minimum number of frames among all input data sets. */
size_t Action_FilterByData::DetermineFrames() const {
  if (Dsets_.empty()) return 0;
  size_t nframes = Dsets_[0]->Size();
  for (Array1D::const_iterator it = Dsets_.begin(); it != Dsets_.end(); ++it)
  {
    if ((*it)->Size() < nframes)
      nframes = (*it)->Size();
  }
  // Warn if any datasets are larger.
  for (Array1D::const_iterator it = Dsets_.begin(); it != Dsets_.end(); ++it)
  {
    if ((*it)->Size() > nframes)
      mprintf("Warning: '%s' size %zu is larger than other sets; only processing %zu\n",
              (*it)->legend(), (*it)->Size(), nframes);
  }
  return nframes;
}

// Action_FilterByData::Print()
void Action_FilterByData::Print() {
  if (!multi_) {
    mprintf("    FILTER: %i frames passed through, %i frames were filtered out.\n",
            Npassed_, Nfiltered_);
    for (unsigned int ds = 0; ds < Dsets_.size(); ds++)
      mprintf("\t%.4f < '%s' < %.4f\n", Min_[ds], Dsets_[ds]->legend(), Max_[ds]);
  }
}
