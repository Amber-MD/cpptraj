#include "Analysis.h"
#include "CpptrajStdio.h"
#include <cstddef> // NULL
// Analysis

// CONSTRUCTOR
Analysis::Analysis() {
  debug = 0;
  analyzeArg = NULL;
  noSetup=false;
}

// DESTRUCTOR
Analysis::~Analysis() {
}

/* Analysis::SetArg()
 * Set analyzeArg to be a copy of the input argument list.
 */
void Analysis::SetArg(ArgList *argIn) {
  analyzeArgs = *argIn;
  analyzeArg = &analyzeArgs;
}

/* Analysis::SetDebug()
 * Set the debug level.
 */
void Analysis::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("ANALYSIS DEBUG LEVEL SET TO %i\n",debug);
}

/* Analysis::AnalysisCommand()
 */
const char *Analysis::AnalysisCommand() {
  return analyzeArg->Command();
}

/* Analysis::CmdLine()
 * Print the command and all args
 */
const char *Analysis::CmdLine() {
  return analyzeArg->ArgLine();
}
