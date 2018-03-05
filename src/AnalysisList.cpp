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

#ifdef MPI
// FIXME Kludge
class Analysis_Placeholder : public Analysis {
  public:
    DispatchObject* Alloc() const { return 0; }
    void Help() const {}
    Analysis_Placeholder() : Analysis(HIDDEN) {}
    RetType Setup(ArgList&, AnalysisSetup&, int) { return Analysis::OK; }
    RetType Analyze() { return Analysis::OK; }
};

/** In parallel currently only trajComm masters will do analysis. This adds
  * a placeholder so that threads can remain in sync.
  * FIXME this is really just a kludge
  */
void AnalysisList::AddPlaceholder() {
  AnaHolder ana;
  ana.ptr_ = new Analysis_Placeholder();
  ana.args_ = ArgList();
  ana.status_ = SETUP;
  analysisList_.push_back( ana );
}
#endif

// AnalysisList::DoAnalyses()
int AnalysisList::DoAnalyses() {
  if (analysisList_.empty()) return 0;
  int err = 0;
  mprintf("\nANALYSIS: Performing %zu analyses:\n",analysisList_.size());
  for (Aarray::const_iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
  {
    if ( ana->status_ == SETUP ) {
      mprintf("  %u: [%s]\n", ana - analysisList_.begin(), ana->args_.ArgLine());
      if (ana->ptr_->Analyze()==Analysis::ERR) {
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
      mprintf("  %u: [%s]\n", ana - analysisList_.begin(), ana->args_.ArgLine());
  }
}
