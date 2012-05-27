// PtrajAnalysis
#include <cstdlib>
#include <cstring>
#include "Analysis_PtrajAnalysis.h"
#include "CpptrajStdio.h"
#include "ptraj_convert.h"

// CONSTRUCTOR
PtrajAnalysis::PtrajAnalysis() {
  argumentStack=NULL;
  analyzeinfo=NULL;
}

// DESTRUCTOR
PtrajAnalysis::~PtrajAnalysis() {
  // Free arguments
  FreeArgumentStack(argumentStack);
  // Free analyzeinfo
  if (analyzeinfo!=NULL) {
    analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_CLEANUP);
    // NOTE: Clean state only for ANALYZE_MODES
    if (analyzeinfo->type==ANALYZE_MODES) 
      FreeState(analyzeinfo->state); 
    free(analyzeinfo);
  }
}

// PtrajAnalysis::Setup()
int PtrajAnalysis::Setup(DataSetList *datasetlist) {
  // Set the global ptraj debug level
  SetPrnlev(debug_);
  return 0;
}

// PtrajAnalysis::Analyze()
/// Call ptraj analyze routine
/** Since ptraj actions are not set up properly until Action::Setup has 
  * been called during trajectory processing, the analyze routines have
  * to be set up here since Analyze::Setup is called before 
  * Action::Setup.
  */
int PtrajAnalysis::Analyze() {
  // Set up analyze information structure
  analyzeinfo = (analyzeInformation*) malloc( sizeof(analyzeInformation));
  INITIALIZE_analyzeInformation( analyzeinfo );
  analyzeinfo->type = ANALYZE_NOOP;

  // Set the analyze type and function based on arg after command
  // (analyze <arg>)
  if        ( analyzeArgs_.ArgIs(1,"correlationcoe")   ) {
    analyzeinfo->type = ANALYZE_CORRELATIONCOEFFICIENT;
    analyzeinfo->fxn  = (analyzeFunction) analyzeCorrelationCoefficient;
  } else if ( analyzeArgs_.ArgIs(1,"crank")            ) {
    analyzeinfo->type = ANALYZE_CRANKSHAFT;
    analyzeinfo->fxn  = (analyzeFunction) analyzeCrankshaft;
  } else if ( analyzeArgs_.ArgIs(1,"stat")             ) {
    analyzeinfo->type = ANALYZE_STATISTICS;
    analyzeinfo->fxn  = (analyzeFunction) analyzeStatistics;
  } else {
    mprinterr("Error: PtrajAnalysis: Unrecognized command: %s\n",analyzeArgs_.ArgAt(1));
  }

  // If still ANALYZE_NOOP there was a problem, free and return
  if (analyzeinfo->type == ANALYZE_NOOP) {
    free( analyzeinfo );
    analyzeinfo = NULL;
    return 1;
  }

  // Mark the arg after 'analyze' command
  analyzeArgs_.MarkArg(1);
  // Convert the remaining args in ArgList to the argStackType used by 
  // functions in ptraj_analyze.c
  argumentStack = CreateArgumentStack(analyzeArgs_, debug_);
  analyzeinfo->carg1 = (void *) &argumentStack;

  // Initialize state memory - only for ANALYZE_MODES
  if (analyzeinfo->type == ANALYZE_MODES) {
    analyzeinfo->state = CreateState( analyzeParm_, analyzeParm_->Nframes() ); // NOTE: maxFrames?
  }

  mprintf("    PTRAJ ANALYZE: [%s %s]\n",analyzeArgs_.Command(),analyzeArgs_.ArgAt(1));

  // NOTE: Eventually need to do something to set up scalarStack

  // Call setup
  if ( analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_SETUP) < 0 )
    return 1;

  // Print status
  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_STATUS);

  // Set prnlev

  // Analyze
  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_ACTION);

  return 0;
}

// PtrajAnalysis::Print()
void PtrajAnalysis::Print(DataFileList *datafilelist) {

  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_PRINT);

}
 
