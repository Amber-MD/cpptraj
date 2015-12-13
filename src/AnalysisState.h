#ifndef INC_ANALYSISSTATE_H
#define INC_ANALYSISSTATE_H
#include "ActionState.h"
/* \file AnalysisState.h
   \brief Classes used to wrap arguments passed to Analyses.
 */
/** The AnalysisSetup class is used to pass in the master DataSetList and
  * DataFileList. Since currently the requirements for Analysis::Setup()
  * are the same as Action::Init(), use the ActionInit class.
  */
typedef ActionInit AnalysisSetup;
#endif
