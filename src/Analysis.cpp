#include "Analysis.h"
#include "CpptrajStdio.h"
#include <cstddef> // NULL
// Analysis

// CONSTRUCTOR
Analysis::Analysis() :
  debug_(0),
  analyzeParm_(NULL),
  isSetup_(false)
{}

// DESTRUCTOR
Analysis::~Analysis() { }

// Analysis::SetArg()
/** Set analyzeArg to be a copy of the input argument list.
  */
void Analysis::SetArg(const ArgList &argIn) {
  analyzeArgs_ = argIn;
}

// Analysis::SetDebug()
/** Set the debug level.  */
void Analysis::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("ANALYSIS DEBUG LEVEL SET TO %i\n",debug_);
}

// Analysis::SetParm()
void Analysis::SetParm(TopologyList *parmfilelist) {
  analyzeParm_ = parmfilelist->GetParm(analyzeArgs_);
}

// Analysis::AnalysisCommand()
const char *Analysis::AnalysisCommand() {
  return analyzeArgs_.Command();
}

// Analysis::CmdLine()
/** Print the command and all args */
const char *Analysis::CmdLine() {
  return analyzeArgs_.ArgLine();
}
