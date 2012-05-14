#include "AnalysisList.h"
#include "CpptrajStdio.h"
// All analysis classes go here
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_PtrajAnalysis.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"

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
 
  if      (argIn.CommandIs("histogram")) { Ana = new Hist(); }
  else if (argIn.CommandIs("hist"))      { Ana = new Hist(); }
  else if (argIn.CommandIs("corr"))      { Ana = new Corr(); }
  else if (argIn.CommandIs("analyze")  ) { 
    if (argIn.ArgAt(1) == NULL) return 1;
    if (argIn[1] == "matrixtest")
      Ana = new Analysis_Matrix;
    else if (argIn[1] == "timecorrtest")
      Ana = new Analysis_Timecorr;
    else
      Ana = new PtrajAnalysis(); 
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
  mprintf("ANALYSIS: Setting up %zu analyses:\n",analysisList_.size());
  int iana = 0;
  for (aListIt ana = analysisList_.begin(); ana != analysisList_.end(); ++ana) {
    mprintf("  %i: [%s]\n",iana++, (*ana)->CmdLine());
    (*ana)->SetSetup(true);
    if ((*ana)->Setup(datasetlist)) {
      mprinterr("Error setting up analysis %i [%s] - skipping.\n",iana,
                (*ana)->AnalysisCommand());
      (*ana)->SetSetup(false);
      ++nfail;
    }
    (*ana)->SetParm(parmfilelist);
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
