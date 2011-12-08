#include "Analysis.h"
#include "CpptrajStdio.h"
#include <cstddef> // NULL
// Analysis

// CONSTRUCTOR
Analysis::Analysis() {
  debug = 0;
  noSetup=false;
  analyzeParm = NULL;
}

// DESTRUCTOR
Analysis::~Analysis() {
}

// Analysis::SetArg()
/** Set analyzeArg to be a copy of the input argument list.
  */
void Analysis::SetArg(const ArgList &argIn) {
  analyzeArgs = argIn;
}

// Analysis::SetDebug()
/** Set the debug level.  */
void Analysis::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("ANALYSIS DEBUG LEVEL SET TO %i\n",debug);
}

// Analysis::SetParm()
void Analysis::SetParm(ParmFileList *parmfilelist) {
  analyzeParm = parmfilelist->GetParm(analyzeArgs);
}

// Analysis::AnalysisCommand()
const char *Analysis::AnalysisCommand() {
  return analyzeArgs.Command();
}

// Analysis::CmdLine()
/** Print the command and all args */
const char *Analysis::CmdLine() {
  return analyzeArgs.ArgLine();
}
