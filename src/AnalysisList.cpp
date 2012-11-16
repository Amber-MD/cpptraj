#include "AnalysisList.h"
#include "CpptrajStdio.h"
// All analysis classes go here
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"
#include "Analysis_IRED.h"
#include "Analysis_Modes.h"
#include "Analysis_CrankShaft.h"
#include "Analysis_Statistics.h"
#include "Analysis_CrossCorr.h"
#include "Analysis_AutoCorr.h"
#include "Analysis_Lifetime.h"
#include "Analysis_FFT.h"
#include "Analysis_CrdFluct.h"

const DispatchObject::Token AnalysisList::DispatchArray[] = {
  { DispatchObject::ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, 0 },
  { DispatchObject::ANALYSIS, "crdfluct", Analysis_CrdFluct::Alloc, Analysis_CrdFluct::Help, 0 },
  { DispatchObject::ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "diagmatrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, 0 },
  { DispatchObject::ANALYSIS, "fft", Analysis_FFT::Alloc, Analysis_FFT::Help, 0 },
  { DispatchObject::ANALYSIS, "hist", Analysis_Hist::Alloc, Analysis_Hist::Help, 0 },
  { DispatchObject::ANALYSIS, "histogram", Analysis_Hist::Alloc, Analysis_Hist::Help, 0 },
  { DispatchObject::ANALYSIS, "ired", Analysis_IRED::Alloc, Analysis_IRED::Help, 0 },
  { DispatchObject::ANALYSIS, "lifetime", Analysis_Lifetime::Alloc, Analysis_Lifetime::Help, 0 },
  { DispatchObject::ANALYSIS, "matrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, 0 },
  { DispatchObject::ANALYSIS, "modes", Analysis_Modes::Alloc, Analysis_Modes::Help, 0 },
  { DispatchObject::ANALYSIS, "timecorr", Analysis_Timecorr::Alloc, Analysis_Timecorr::Help, 0 },
  { DispatchObject::ANALYSIS, "stat", Analysis_Statistics::Alloc, Analysis_Statistics::Help, 0 },
  { DispatchObject::NONE,        0,                  0,                 0, 0 }
};

// CONSTRUCTOR
AnalysisList::AnalysisList() :
  debug_(0)
{}

// DESTRUCTOR
AnalysisList::~AnalysisList() {
  Clear();
}

void AnalysisList::Clear() {
  for (aListType::iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
    delete *ana;
  analysisList_.clear();
}

// AnalysisList::SetDebug()
/** Set Analysis list debug level. */
void AnalysisList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("AnalysisList DEBUG LEVEL SET TO %i\n",debug_);
}

// AnalysisList::AddAnalysis()
/** Add specified analysis to the analysis list with given args and 
  * DataSetList.
  */
int AnalysisList::AddAnalysis(DispatchObject::DispatchAllocatorType Alloc, ArgList& argIn,
                              TopologyList* PFLin, DataSetList* DSLin)
{
  Analysis* ana = (Analysis*)Alloc();
  // Attempt to set up analysis
  if (ana->Setup( argIn, DSLin, PFLin, debug_) != Analysis::OK) {
    mprinterr("Error: Could not setup analysis [%s]\n", argIn.Command());
    return 1;
  }
  argIn.CheckForMoreArgs();
  analysisList_.push_back( ana );
  analysisCmd_.push_back( argIn.ArgLine() );
  analysisStatus_.push_back( SETUP );
  return 0;
}

// AnalysisList::DoAnalyses()
void AnalysisList::DoAnalyses(DataFileList *datafilelist) {
  if (analysisList_.empty()) return;
  mprintf("\nANALYSIS: Performing %zu analyses:\n",analysisList_.size());
  unsigned int ananum = 0;
  for (aListType::iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana) {
    if ( analysisStatus_[ananum] == SETUP ) {
      mprintf("  %u: [%s]\n", ananum, analysisCmd_[ananum].c_str());
      if ((*ana)->Analyze()==Analysis::OK) 
        (*ana)->Print(datafilelist); 
      // NOTE: Move print function ??
    }
    ++ananum;
  }
  mprintf("\n");
  //mprintf("    ...................................................\n\n");
}

void AnalysisList::List() {
  unsigned int ananum = 0;
  for (aListType::iterator ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
  {
    mprintf("  %u: [%s]\n", ananum, analysisCmd_[ananum].c_str());
    ++ananum;
  }
}
