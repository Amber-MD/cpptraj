#include "DataFilter.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
DataFilter::DataFilter() :
  filterSet_(0),
  Npassed_(0),
  Nfiltered_(0),
  debug_(0),
  multi_(false)
{}

/** Keywords recognized by InitFilter(). */
const char* DataFilter::Keywords() {
  return "";
}

/** Process arguments, get/set up data sets. */
int DataFilter::InitFilter(ArgList& argIn, DataSetList& DSL, DataFileList& DFL, int debugIn) {
  filterSet_ = 0;
  Max_.clear();
  Min_.clear();
  inpSets_.clear();
  outSets_.clear();
  Npassed_ = 0;
  Nfiltered_ = 0;
  debug_ = debugIn;

  multi_ = argIn.hasKey("multi");
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())
    dsname = DSL.GenerateDefaultName("Filter");
  DataFile* maxminfile = DFL.AddDataFile( argIn.GetStringKey("out"), argIn );
  // Get min and max args.
  while (argIn.Contains("min"))
    Min_.push_back( argIn.getKeyDouble("min", 0.0) );
  while (argIn.Contains("max"))
    Max_.push_back( argIn.getKeyDouble("max", 0.0) );
  if (Min_.empty()) {
    mprinterr("Error: At least one 'min' arg must be specified.\n");
    return 1;
  }
  if (Max_.empty()) {
    mprinterr("Error: At least one 'max' arg must be specified.\n");
    return 1;
  }
  if (Min_.size() != Max_.size()) {
    mprinterr("Error: # of 'min' args (%zu) != # of 'max' args (%zu)\n",
              Min_.size(), Max_.size());
    return 1;
  }
  // Get DataSets from remaining arguments
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    DataSetList sets = DSL.GetMultipleSets( dsarg );
    if (sets.empty())
      mprintf("Warning: '%s' selects no sets.\n", dsarg.c_str());
    else {
      for (DataSetList::const_iterator ds = sets.begin(); ds != sets.end(); ++ds)
        inpSets_.push_back( *ds );
    }
  }
  if (inpSets_.empty()) {
    mprinterr("Error: No data sets specified.\n");
    return 1;
  }

  // Check number of DataSets and min/max args.
  if ( inpSets_.size() < Min_.size() ) {
    mprinterr("Error: More 'min'/'max' args (%zu) than data sets (%zu).\n",
              Min_.size(), inpSets_.size());
    return 1;
  }
  if ( inpSets_.size() > Min_.size() ) {
    unsigned int Nremaining = inpSets_.size() - Min_.size();
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
    filterSet_ = (DataSet_integer*)DSL.AddSet( DataSet::INTEGER, dsname );
    if (filterSet_ == 0) return 1;
    if (maxminfile != 0)
      maxminfile->AddDataSet( (DataSet*)filterSet_ );
  } else {
    for (unsigned int idx = 0; idx < inpSets_.size(); idx++) {
      DataSet* ds = DSL.AddSet(DataSet::INTEGER, MetaData(dsname, idx));
      if (ds == 0) return 1;
      ds->SetLegend("Filter("+ds->Meta().PrintName()+")");
      outSets_.push_back( ds );
      if (maxminfile != 0)
        maxminfile->AddDataSet( ds );
    }
  }

  return 0;
}
