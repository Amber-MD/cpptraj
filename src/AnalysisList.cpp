#include "AnalysisList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AnalysisList::AnalysisList() : debug_(0) {}

// DESTRUCTOR
AnalysisList::~AnalysisList() {
  Clear();
}

void AnalysisList::Clear() {
  for (aListType::iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
    delete *ana;
  analysisList_.clear();
  analysisCmd_.clear();
  analysisStatus_.clear();
}

// AnalysisList::SetDebug()
/** Set Analysis list debug level. */
void AnalysisList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("AnalysisList DEBUG LEVEL SET TO %i\n",debug_);
}

// AnalysisList::DoAnalyses()
void AnalysisList::DoAnalyses() {
  if (analysisList_.empty()) return;
  mprintf("\nANALYSIS: Performing %zu analyses:\n",analysisList_.size());
  unsigned int ananum = 0;
  for (aListType::iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana) {
    if ( analysisStatus_[ananum] == SETUP ) {
      mprintf("  %u: [%s]\n", ananum, analysisCmd_[ananum].c_str());
      if ((*ana)->Analyze()==Analysis::ERR)
        mprinterr("Error: in Analysis # %u\n", ananum); 
    }
    ++ananum;
  }
  mprintf("\n");
  //mprintf("    ...................................................\n\n");
}

void AnalysisList::List() const {
  mprintf("ANALYSES:\n");
  if (analysisList_.empty())
    mprintf("  No Analyses.\n");
  else
    for (unsigned int ananum = 0; ananum < analysisList_.size(); ++ananum)
      mprintf("  %u: [%s]\n", ananum, analysisCmd_[ananum].c_str());
}
