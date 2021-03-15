#include "DataFilter.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_2D.h"
#include "DataSet_3D.h"
#include "DataSet_integer.h"
#include "DataFile.h"

/** CONSTRUCTOR */
DataFilter::DataFilter() :
  filterSet_(0),
  SetToBeFiltered_(0),
  FilteredSet_(0),
  npassedSet_(0),
  nfilteredSet_(0),
  maxminfile_(0),
  countFile_(0),
  masterDSL_(0),
  debug_(0),
  multi_(false)
{
  Nresult_[PASSED] = 0;
  Nresult_[FILTERED] = 0;
}

/** Keywords recognized by InitFilter(). */
void DataFilter::PrintKeywords() {
  mprintf("\t{<dataset arg> min <min> max <max> ...} [out <file>] [name <setname>]\n"
          "\t{[multi] | [filterset <set> [newset <newname>]] [countout <countfile>]}\n");
}

/** A more detailed description of some keywords. */
void DataFilter::PrintKeywordDescriptions() {
  mprintf("  If 'multi' is specified then only filter data sets will be created for each\n"
          "  data set instead.\n"
          "  If 'filterset' is specified, the specified <set> will be modified\n"
          "  to only contain '1' frames; cannot be used with 'multi'. If 'newset'\n"
          "  is also specified, a new set will be created containing the '1' frames instead.\n"
          "  The 'filterset' functionality only works for 1D scalar sets.\n"
          "  If 'countout' is specified, the final number of elements passed and filtered\n"
          "  out will be written to <countfile>.\n");
}

/** Process arguments, get/set up data sets. */
int DataFilter::InitFilter(ArgList& argIn, DataSetList& DSL, DataFileList& DFL, int debugIn) {
  filterSet_ = 0;
  SetToBeFiltered_ = 0;
  FilteredSet_ = 0;
  npassedSet_ = 0;
  nfilteredSet_ = 0;
  maxminfile_ = 0;
  countFile_ = 0;
  Xvals_.clear();
  Max_.clear();
  Min_.clear();
  inpSets_.clear();
  outSets_.clear();
  Nresult_[PASSED] = 0;
  Nresult_[FILTERED] = 0;
  debug_ = debugIn;
  outIdx_ = 0;

  masterDSL_ = &DSL;

  multi_ = argIn.hasKey("multi");

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())
    dsname = DSL.GenerateDefaultName("Filter");
  maxminfile_ = DFL.AddDataFile( argIn.GetStringKey("out"), argIn );
  countFile_ = DFL.AddDataFile( argIn.GetStringKey("countout"), argIn );
  // Get set to be filtered
  std::string nameOfSetToBeFiltered = argIn.GetStringKey("filterset");
  std::string resultingSetName = argIn.GetStringKey("newset");

  // Check for options incompatible with 'multi'.
  if (multi_) {
    if (!nameOfSetToBeFiltered.empty()) {
      mprinterr("Error: 'filterset' can not be used with 'multi' keyword.\n");
      return 1;
    }
    if (countFile_ != 0) {
      mprinterr("Error: 'countout' can not be used with 'multi' keyword.\n");
      return 1;
    }
  }
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
  // Get input DataSets from remaining arguments
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    DataSetList sets = DSL.GetMultipleSets( dsarg );
    if (sets.empty())
      mprintf("Warning: '%s' selects no sets.\n", dsarg.c_str());
    else {
      for (DataSetList::const_iterator ds = sets.begin(); ds != sets.end(); ++ds)
      {
        if ((*ds)->Group() != DataSet::SCALAR_1D &&
            (*ds)->Group() != DataSet::MATRIX_2D &&
            (*ds)->Group() != DataSet::GRID_3D)
        {
          mprinterr("Error: '%s' is not a 1D, 2D, or 3D data set.", (*ds)->legend());
          return 1;
        }
        inpSets_.push_back( *ds );
      }
    }
    dsarg = argIn.GetStringNext();
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
    if (maxminfile_ != 0)
      maxminfile_->AddDataSet( (DataSet*)filterSet_ );
    // Create count sets
    MetaData count_md(dsname, "npassed");
    npassedSet_ = DSL.AddSet( DataSet::UNSIGNED_INTEGER, count_md );
    count_md.SetAspect("nfiltered");
    nfilteredSet_ = DSL.AddSet( DataSet::UNSIGNED_INTEGER, count_md );
    if (npassedSet_ == 0 || nfilteredSet_ == 0) {
      mprinterr("Error: Could not set up count sets.\n");
      return 1;
    }
    if (countFile_ != 0) {
      countFile_->AddDataSet( npassedSet_ );
      countFile_->AddDataSet( nfilteredSet_ );
    }
    // Get set to be filtered if specified.
    if (!nameOfSetToBeFiltered.empty()) {
      DataSet* ds = DSL.GetDataSet( nameOfSetToBeFiltered );
      if (ds == 0) {
        mprinterr("Error: Set to be filtered '%s' not found.\n", nameOfSetToBeFiltered.c_str());
        return 1;
      }
      if (ds->Group() != DataSet::SCALAR_1D)
      {
        mprinterr("Error: '%s' is not 1D scalar.\n", ds->legend());
        return 1;
      }
      //if (ds->Size() < nframes) { TODO re-enable this check in e.g. Setup
      //  mprinterr("Error: Set to be filtered size (%zu) is too small (%zu)\n",
      //            ds->Size(), nframes);
      //  return 1;
      //}
      SetToBeFiltered_ = (DataSet_1D*)ds;
      // Create new filter set
      FilteredSet_ = (DataSet_1D*)DataSetList::Allocate(SetToBeFiltered_->Type());
      MetaData md(SetToBeFiltered_->Meta());
      if (!resultingSetName.empty()) {
        md.SetName(resultingSetName);
        ds = DSL.CheckForSet( md );
        if (ds != 0) {
          mprinterr("Error: New set name '%s' already exists.\n", resultingSetName.c_str());
          return 1;
        }
        mprintf("\tA new filtered set will be created from set '%s'\n", SetToBeFiltered_->legend());
      } else
        mprintf("\tFiltering set '%s'\n", SetToBeFiltered_->legend());
      FilteredSet_->SetMeta( md );
      FilteredSet_->SetDim(0, SetToBeFiltered_->Dim(0));
      // Do not allocate here to save memory.
    }
  } else {
    for (unsigned int idx = 0; idx < inpSets_.size(); idx++) {
      DataSet* ds = DSL.AddSet(DataSet::INTEGER, MetaData(dsname, idx));
      if (ds == 0) return 1;
      ds->SetLegend("Filter("+ds->Meta().PrintName()+")");
      outSets_.push_back( ds );
      if (maxminfile_ != 0)
        maxminfile_->AddDataSet( ds );
    }
  }

  return 0;
}

/** \return Minimum number of elements among all input data sets. */
size_t DataFilter::MinNumElements() const {
  if (inpSets_.empty()) return 0;
  size_t nelements = inpSets_[0]->Size();
  for (SetArray::const_iterator it = inpSets_.begin(); it != inpSets_.end(); ++it)
  {
    if ((*it)->Size() < nelements)
      nelements = (*it)->Size();
  }
  // Warn if any datasets are larger.
  for (SetArray::const_iterator it = inpSets_.begin(); it != inpSets_.end(); ++it)
  {
    if ((*it)->Size() > nelements)
      mprintf("Warning: '%s' size %zu is larger than other sets; only processing %zu\n",
              (*it)->legend(), (*it)->Size(), nelements);
  }
  return nelements;
}


/** Get input value for specified set index. */
double DataFilter::GetInpValue(unsigned int setIdx, unsigned int inpIdx) const {
  double dVal = 0;
  if (inpSets_[setIdx]->Group() == DataSet::SCALAR_1D) {
    dVal = ((DataSet_1D*)inpSets_[setIdx])->Dval(inpIdx);
  } else if (inpSets_[setIdx]->Group() == DataSet::MATRIX_2D) {
    dVal = ((DataSet_2D*)inpSets_[setIdx])->GetElement(inpIdx);
  } else if (inpSets_[setIdx]->Group() == DataSet::GRID_3D) {
    dVal = (*((DataSet_3D*)inpSets_[setIdx]))[inpIdx];
  } else {
    mprinterr("Error: Unhandled set type in DataFilter::GetInpValue().\n");
  }
  return dVal;
}

/** Integer values to be put into result data set(s). */
const int DataFilter::ResultValue[2] = {
  1, // PASSED
  0, // FILTERED
};

/** Filter the specified index.
  * \return 1 if index was filtered, 0 if index passed.
  */
DataFilter::ResultType DataFilter::FilterIndex(unsigned int inpIdx) {
  ResultType result = PASSED;
  if (!multi_) {
    // Check if frame is within max/min of every data set.
    for (unsigned int idx = 0; idx != inpSets_.size(); ++idx)
    {
      double dVal = GetInpValue(idx, inpIdx);
      //mprintf("DBG: maxmin[%u]: dVal = %f, min = %f, max = %f\n",idx,dVal,Min_[idx],Max_[idx]);
      // If value from dataset not within min/max, exit now.
      if (dVal < Min_[idx] || dVal > Max_[idx]) {
        result = FILTERED; 
        break;
      }
    }
    filterSet_->Add( outIdx_, ResultValue + (int)result );
    // If filtering set and element passed, add element to new set (FilteredSet).
    if (SetToBeFiltered_ != 0 && result == PASSED) {
      Xvals_.push_back( SetToBeFiltered_->Xcrd(inpIdx));
      FilteredSet_->Add(Nresult_[PASSED], SetToBeFiltered_->VoidPtr(inpIdx));
    }
    Nresult_[result]++;
  } else {
    // For each data set 1 if within min/max, 0 otherwise.
    for (unsigned int idx = 0; idx != inpSets_.size(); ++idx)
    {
      double dVal = GetInpValue(idx, inpIdx);
      if (dVal < Min_[idx] || dVal > Max_[idx]) 
        outSets_[idx]->Add( outIdx_, ResultValue + (int)FILTERED );
      else
        outSets_[idx]->Add( outIdx_, ResultValue + (int)PASSED );
    }
  }
  outIdx_++;
  return result;
}

/** Print input sets and min/max values to stdout. */
void DataFilter::PrintInputSets() const {
  for (unsigned int idx = 0; idx < inpSets_.size(); idx++)
    mprintf("\t%.4f < '%s' < %.4f\n", Min_[idx], inpSets_[idx]->legend(), Max_[idx]);
}

/** Perform any actions necessary to finish filtering. */
int DataFilter::Finalize() const {
  if (!multi_) {
    unsigned int uval = Npassed();
    npassedSet_->Add(0, &uval);
    uval = Nfiltered();
    nfilteredSet_->Add(0, &uval);
  }

  if (SetToBeFiltered_ != 0) {
    if (SetToBeFiltered_->Meta().Match_Exact(FilteredSet_->Meta()))
      masterDSL_->RemoveSet( SetToBeFiltered_ );
    DataSetList::DataListType inputSets(1);
    inputSets[0] = FilteredSet_;
    masterDSL_->AddOrAppendSets( FilteredSet_->Dim(0).Label(), Xvals_, inputSets );
    // AddOrAppendSets() will free memory.
    //FilteredSet_ = 0;
  }
  return 0;
}
