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

const DispatchObject::Token AnalysisList::DispatchArray[] = {
  { DispatchObject::ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, 0 },
  { DispatchObject::ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, 0 },
  { DispatchObject::ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, 0 },
  { DispatchObject::ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, 0 },
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
  for (aListIt ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
    delete *ana;
}

// AnalysisList::SetDebug()
/** Set Analysis list debug level. */
void AnalysisList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("AnalysisList DEBUG LEVEL SET TO %i\n",debug_);
}

int AnalysisList::AddAnalysis(DispatchObject::DispatchAllocatorType Alloc, ArgList const& argIn)
{
  Analysis* ana = (Analysis*)Alloc();
  ana->SetArg( argIn );
  ana->SetDebug( debug_ );
  analysisList_.push_back( ana );
  return 0;
}

// AnalysisList::Setup()
/** Set up all analysis in list with given datasetlist. Also set the parm
  * (first parm will be set if parm/parmindex keywords not specified).
  */
int AnalysisList::Setup(DataSetList *datasetlist, TopologyList *parmfilelist) {
  int nfail = 0;
  if (analysisList_.empty()) return 0;
  mprintf("\nANALYSIS: Setting up %zu analyses:\n",analysisList_.size());
  int iana = 0;
  for (aListIt ana = analysisList_.begin(); ana != analysisList_.end(); ++ana) {
    // Set parm for analysis.
    (*ana)->SetParm(parmfilelist);
    mprintf("  %i: [%s] (Parm: %s)\n", iana++, (*ana)->CmdLine(), (*ana)->ParmName());
    (*ana)->SetSetup(true);
    if ((*ana)->Setup(datasetlist)) {
      mprinterr("Error setting up analysis %i [%s] - skipping.\n",iana,
                (*ana)->AnalysisCommand());
      (*ana)->SetSetup(false);
      ++nfail;
    }
    (*ana)->CheckForMoreArgs();
  }
  mprintf("\n");   
  //mprintf("    ...................................................\n\n");
  return nfail;
}

// AnalysisList::Analyze()
void AnalysisList::Analyze(DataFileList *datafilelist) {
  if (analysisList_.empty()) return;
  mprintf("\nANALYSIS: Performing %zu analyses:\n",analysisList_.size());
  int iana = 0;
  for (aListIt ana = analysisList_.begin(); ana != analysisList_.end(); ++ana) {
    if ((*ana)->IsSetup()) {
      mprintf("  %i: [%s]\n",iana, (*ana)->CmdLine());
      if ((*ana)->Analyze()==0) 
        (*ana)->Print(datafilelist); 
      // NOTE: Move print function ??
    }
    ++iana;
  }
  mprintf("\n");
  //mprintf("    ...................................................\n\n");
}

void AnalysisList::List() {
  unsigned int iana = 0;
  for (aListIt ana = analysisList_.begin(); ana != analysisList_.end(); ++ana)
    mprintf("  %u: [%s] (Parm: %s)\n", iana++, (*ana)->CmdLine(), (*ana)->ParmName());
}
