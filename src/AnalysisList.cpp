#include "AnalysisList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AnalysisList::AnalysisList() : debug_(0) {}

// DESTRUCTOR
AnalysisList::~AnalysisList() { Clear(); }

// AnalysisList::Clear()
void AnalysisList::Clear() {
  for (Aarray::const_iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
    delete ana->ptr_;
  analysisList_.clear();
}

// AnalysisList::AddAnalysis()
/** Add specified analysis to the analysis list with given args and 
  * DataSetList.
  */
int AnalysisList::AddAnalysis(Analysis* anaIn, ArgList& argIn, AnalysisSetup& setup)
{
  if (anaIn == 0) {
    mprinterr("Internal Error: AddAnalysis() called with null Analysis.\n");
    return 1;
  }
  AnaHolder ana;
  ana.ptr_ = anaIn; 
  ana.args_ = argIn;
  // Attempt to set up analysis
  if (ana.ptr_->Setup( argIn, setup, debug_) != Analysis::OK) {
    mprinterr("Error: Could not setup analysis [%s]\n", argIn.Command());
    delete ana.ptr_;
    return 1;
  }
  ana.status_ = SETUP;
  analysisList_.push_back( ana );
  if (argIn.CheckForMoreArgs()) return 1;
  return 0;
}

// AnalysisList::DoAnalyses()
int AnalysisList::DoAnalyses() {
  if (analysisList_.empty()) return 0;
  int err = 0;
  mprintf("\nANALYSIS: Performing %zu analyses:\n",analysisList_.size());
  for (Aarray::const_iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
  {
    if ( ana->status_ == SETUP ) {
      mprintf("  %li: [%s]\n", ana - analysisList_.begin(), ana->args_.ArgLine());
      Analysis::RetType ret;
#     ifdef MPI
      if (ana->ptr_->IsParallel()) {
        ret = ana->ptr_->Analyze();
      } else {
        mprintf("Warning: Analysis '%s' does not currently use multiple MPI processes.\n", ana->args_.Command());
        if (Parallel::TrajComm().Master())
          ret = ana->ptr_->Analyze();
        else
          ret = Analysis::OK;
        Parallel::World().MasterBcast( &ret, 1, MPI_INT );
      }
      int err;
      if (ret == Analysis::ERR) {
        rprinterr("Error: In parallel, analysis '%s' failed.\n", ana->args_.Command());
        err = 1;
      } else
        err = 0;
      if (Parallel::World().CheckError( err ))
        ret = Analysis::ERR;
#     else /* MPI */
      ret = ana->ptr_->Analyze();
#     endif /* MPI */
      if (ret == Analysis::ERR) {
        mprinterr("Error: In Analysis [%s]\n", ana->args_.Command()); // TODO exit? Set INACTIVE?
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

// AnalysisList::List()
void AnalysisList::List() const {
  if (!analysisList_.empty()) {
    mprintf("\nANALYSES (%zu total):\n", analysisList_.size());
    for (Aarray::const_iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
      mprintf("  %li: [%s]\n", ana - analysisList_.begin(), ana->args_.ArgLine());
  }
}
