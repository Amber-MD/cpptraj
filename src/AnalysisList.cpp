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

// AnalysisList::AddAnalysis()
/** Add specific type of analysis to the list.  */
int AnalysisList::AddAnalysis(ArgList &argIn) {
  Analysis *Ana=NULL;
 
  if      (argIn.CommandIs("histogram")) { Ana = new Analysis_Hist(); }
  else if (argIn.CommandIs("hist"))      { Ana = new Analysis_Hist(); }
  else if (argIn.CommandIs("corr"))      { Ana = new Analysis_Corr(); }
  else if (argIn.CommandIs("crosscorr")) { Ana = new Analysis_CrossCorr(); }
  else if (argIn.CommandIs("autocorr"))  { Ana = new Analysis_AutoCorr(); }
  else if (argIn.CommandIs("lifetime"))  { Ana = new Analysis_Lifetime(); }
  else if (argIn.CommandIs("fft"))       { Ana = new Analysis_FFT(); }
  else if (argIn.CommandIs("analyze")  ) { 
    if (argIn.ArgAt(1) == NULL) return 1;
    if (argIn[1] == "matrix")
      Ana = new Analysis_Matrix;
    else if (argIn[1] == "timecorr")
      Ana = new Analysis_Timecorr;
    else if (argIn[1] == "ired")
      Ana = new Analysis_IRED;
    else if (argIn[1] == "modes")
      Ana = new Analysis_Modes;
    else if (argIn[1] == "crank")
      Ana = new Analysis_CrankShaft;
    else if (argIn[1] == "correlationcoe")
      Ana = new Analysis_Corr; // For backwards compatibility with PTRAJ
    else if (argIn[1] == "stat")
      Ana = new Analysis_Statistics;
    else
      return 1; 
  }
  else return 1;

  // Pass in the argument list
  Ana->SetArg(argIn);

  // Set the debug level
  if (debug_>0) Ana->SetDebug(debug_);

  // Store in analysis list
  analysisList_.push_back(Ana);

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
