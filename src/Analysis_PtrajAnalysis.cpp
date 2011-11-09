// PtrajAnalysis
#include <cstdlib>
#include <cstring>
#include "Analysis_PtrajAnalysis.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
PtrajAnalysis::PtrajAnalysis() {
  argumentStack=NULL;
  analyzeinfo=NULL;
}

// DESTRUCTOR
PtrajAnalysis::~PtrajAnalysis() {
  // Free arguments
  if (argumentStack!=NULL) {
    // free individual args?
    for (int arg=0; arg < argumentStack->nargs; arg++)
      free(argumentStack->arglist[arg]);
    free(argumentStack->arglist);
    free(argumentStack->marked);
    free(argumentStack);
  }
  // Free analyzeinfo
  if (analyzeinfo!=NULL) {
    // NOTE: Only for ANALYZE_MODES
    if (analyzeinfo->type==ANALYZE_MODES) 
      free(analyzeinfo->state); 
    free(analyzeinfo);
  }
}

// PtrajAnalysis::Setup()
int PtrajAnalysis::Setup(DataSetList *datasetlist) {
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
  } else if ( analyzeArgs.ArgIs(1,"matrix")           ) {
    analyzeinfo->type = ANALYZE_MATRIX;
    analyzeinfo->fxn  = (analyzeFunction) analyzeMatrix;
  } else {
    mprinterr("Error: PtrajAnalysis: Unrecognized command: %s\n",analyzeArgs.ArgAt(1));
    return 1;
  }

  // Mark the arg after command
  analyzeArgs.MarkArg(1);
  // Convert the ArgList to the argStackType used by functions in ptraj_analyze.c
  argumentStack = (argStackType*) malloc( sizeof(argStackType) );
  //mprintf("DEBUG:\targumentStack address (init): %x\n",argumentStack);
  argumentStack->arglist = NULL;
  argumentStack->marked = NULL;
  int nargs = 0;
  char *currentArg = analyzeArgs.getNextString();
  while (currentArg != NULL) {
    argumentStack->arglist = (char**) realloc(argumentStack->arglist, (nargs+1) * sizeof(char*));
    argumentStack->arglist[nargs] = (char*) malloc( (strlen(currentArg)+1) * sizeof(char) );
    strcpy(argumentStack->arglist[nargs], currentArg);
    currentArg = analyzeArgs.getNextString();
    nargs++;
  } 
  argumentStack->nargs = nargs;
  argumentStack->marked = (char*) malloc(nargs * sizeof(char));
  memset(argumentStack->marked, 'F', nargs);
  analyzeinfo->carg1 = (void *) &argumentStack;
  //mprintf("DEBUG:\taction->carg1 address (init): %x\n",actioninfo->carg1);
  // NOTE: Should ptraj be freeing up the args?
  // DEBUG
  //argStackType **argumentStackPointer = (argStackType **) actioninfo->carg1;
  //mprintf("DEBUG:\targumentStack address (init 2): %x\n",*argumentStackPointer);
  if (debug>0) printArgumentStack(&argumentStack);

  // Initialize state memory - only for ANALYZE_MODES
  if (analyzeinfo->type == ANALYZE_MODES) {
    analyzeinfo->state = (ptrajState*) malloc( sizeof(ptrajState) );
    INITIALIZE_ptrajState( analyzeinfo->state );
  }

  mprintf("    PTRAJ ANALYZE: [%s %s]\n",analyzeArgs.Command(),analyzeArgs.ArgAt(1));

  // NOTE: Eventually need to do something to set up scalarStack

  // Call setup
  if ( analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_SETUP) < 0 )
    return 1;
  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_STATUS);

  // Set prnlev

  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_ACTION);

  return 0;
}

void PtrajAnalysis::Print(DataFileList *datafilelist) {

  analyzeinfo->fxn(analyzeinfo, NULL, PTRAJ_PRINT);

}
 
