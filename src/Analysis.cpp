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
  if (analyzeArg!=NULL) delete analyzeArg;
}

/* Analysis::SetArg()
 * Set analyzeArg to be a copy of the input argument list.
 */
void Analysis::SetArg(ArgList *argIn) {
  analyzeArg = argIn->Copy();
}

/* Analysis::SetDebug()
 * Set the debug level.
 */
void Analysis::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("ANALYSIS DEBUG LEVEL SET TO %i\n",debug);
}

/* Analysis::Name()
 */
char *Analysis::Name() {
  return analyzeArg->Command();
}

/* Analysis::CmdLine()
 * Print the command and all args
 */
char *Analysis::CmdLine() {
  return analyzeArg->ArgLine();
}
