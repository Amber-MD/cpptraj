// Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Outtraj::Outtraj() {
  //fprintf(stderr,"Outtraj Con\n");
  min=0.0;
  max=0.0;
  Dset=NULL;
} 

// Outtraj::init()
/** Expected call: outtraj <filename> [ trajout args ] 
  *                        [maxmin <dataset> min <min> max <max>
  */
int Outtraj::init() {
  char *datasetName;
  AmberParm *tempParm;

#ifdef MPI
  mprintf("ERROR: OUTTRAJ currently not functional with MPI.\n");
  return 1;
#endif

  outtraj.SetDebug(debug);
  tempParm = PFL->GetParm(actionArgs);
  if (tempParm==NULL) {
    mprinterr("Error: OUTTRAJ: Could not get parm for %s\n",actionArgs.ArgAt(1));
    return 1;
  }
  if ( outtraj.SetupWrite(NULL,&actionArgs,tempParm,AMBERTRAJ) ) return 1;
  mprintf("    OUTTRAJ:");
  outtraj.PrintInfo(1);
  // If maxmin, get the name of the dataset as well as the max and min values.
  datasetName = actionArgs.getKeyString("maxmin",NULL);
  if (datasetName!=NULL) {
    Dset = DSL->Get(datasetName);
    if (Dset==NULL) {
      mprintf("Error: Outtraj maxmin: Could not get dataset %s\n",datasetName);
      return 1;
    } else {
      // Currently dont allow for string datasets
      if (Dset->Type()==STRING) {
        mprintf("Error: Outtraj maxmin: String dataset (%s) not supported.\n",datasetName);
        return 1;
      }
      max = actionArgs.getKeyDouble("max",0.0);
      min = actionArgs.getKeyDouble("min",0.0);
      mprintf("             maxmin: Printing trajectory frames based on %lf <= %s <= %lf\n",
              min, datasetName, max);
    }
  }

  return 0;
} 

// Outtraj::setup()
// Unneeded, Output trajectory setup done right before first write
//int Outtraj::setup() {
//  return 0;  
//}

// Outtraj::action()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
int Outtraj::action() {
  double dVal;
  int iVal;

  // If dataset defined, check if frame is within max/min
  if (Dset!=NULL) {
    if (Dset->Type() == DOUBLE) {
      if (Dset->Get(&dVal, frameNum)) return 1;
    } else if (Dset->Type() == INT) {
      if (Dset->Get(&iVal, frameNum)) return 1;
      dVal = (double) iVal;
    } else
      return 1;
    //mprintf("DBG: maxmin: dVal = %lf\n",dVal);
    // If value from dataset not within min/max, exit now.
    if (dVal < min || dVal > max) return 0;
  }
  if ( outtraj.WriteFrame(frameNum, currentParm, *currentFrame) != 0 ) 
    return 1;
  return 0;
} 

// Outtraj::print()
/** Close trajectory. Indicate how many frames were actually written.
  */
void Outtraj::print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",actionArgs.ArgAt(1),outtraj.NumFramesProcessed());
  outtraj.EndTraj();
}

