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
  // Clear the stacks. Set to NULL since the stacks are global. This 
  // prevents another PtrajAnalysis from trying to free them again.
  //if (vectorStack!=NULL) clearStack( &vectorStack );
  //if (matrixStack!=NULL) clearStack( &matrixStack );
  //if (modesStack!=NULL) clearStack( &modesStack );
  //vectorStack = NULL;
  //matrixStack = NULL;
  //modesStack = NULL;
}

// PtrajAnalysis::Setup()
int PtrajAnalysis::Setup(DataSetList *datasetlist) {
  // Set the global ptraj debug level
  SetPrnlev(debug);
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
  if        ( analyzeArgs.ArgIs(1,"correlationcoe")   ) {
    analyzeinfo->type = ANALYZE_CORRELATIONCOEFFICIENT;
    analyzeinfo->fxn  = (analyzeFunction) analyzeCorrelationCoefficient;
  } else if ( analyzeArgs.ArgIs(1,"crank")            ) {
    analyzeinfo->type = ANALYZE_CRANKSHAFT;
    analyzeinfo->fxn  = (analyzeFunction) analyzeCrankshaft;
#ifndef NO_PTRAJ_ANALYZE
  } else if ( analyzeArgs.ArgIs(1,"matrix")           ) {
    analyzeinfo->type = ANALYZE_MATRIX;
    analyzeinfo->fxn  = (analyzeFunction) analyzeMatrix;
  } else if ( analyzeArgs.ArgIs(1,"timecorr")         ) {
    analyzeinfo->type = ANALYZE_TIMECORR;
    analyzeinfo->fxn  = (analyzeFunction) analyzeTimecorr;
#else
  } else if ( analyzeArgs.ArgIs(1,"matrix")           ) {
    mprinterr("Error: cpptraj compiled with -DNO_PTRAJ_ANALYZE, 'analyze matrix' unsupported.\n");
  } else if ( analyzeArgs.ArgIs(1,"timecorr")         ) {
    mprinterr("Error: cpptraj compiled with -DNO_PTRAJ_ANALYZE, 'analyze timecorr' unsupported.\n");
#endif
  } else if ( analyzeArgs.ArgIs(1,"modes")            ) {
    analyzeinfo->type = ANALYZE_MODES;
    analyzeinfo->fxn  = (analyzeFunction) analyzeModes;
  } else if ( analyzeArgs.ArgIs(1,"set")              ) {
    analyzeinfo->type = ANALYZE_SET;
    analyzeinfo->fxn  = (analyzeFunction) analyzeSet;
  } else if ( analyzeArgs.ArgIs(1,"stat")             ) {
    analyzeinfo->type = ANALYZE_STATISTICS;
    analyzeinfo->fxn  = (analyzeFunction) analyzeStatistics;
  } else if ( analyzeArgs.ArgIs(1,"test")             ) {
    analyzeinfo->type = ANALYZE_TEST;
    analyzeinfo->fxn  = (analyzeFunction) analyzeTest;
  } else {
    mprinterr("Error: PtrajAnalysis: Unrecognized command: %s\n",analyzeArgs.ArgAt(1));
  }

  // If still ANALYZE_NOOP there was a problem, free and return
  if (analyzeinfo->type == ANALYZE_NOOP) {
    free( analyzeinfo );
    analyzeinfo = NULL;
    return 1;
  }

  // Mark the arg after 'analyze' command
  analyzeArgs.MarkArg(1);
  // Convert the remaining args in ArgList to the argStackType used by 
  // functions in ptraj_analyze.c
  argumentStack = CreateArgumentStack(analyzeArgs, debug);
  analyzeinfo->carg1 = (void *) &argumentStack;

  // Initialize state memory - only for ANALYZE_MODES
  if (analyzeinfo->type == ANALYZE_MODES) {
    analyzeinfo->state = CreateState( analyzeParm, analyzeParm->parmFrames ); // NOTE: maxFrames?
  }

  mprintf("    PTRAJ ANALYZE: [%s %s]\n",analyzeArgs.Command(),analyzeArgs.ArgAt(1));

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
 
