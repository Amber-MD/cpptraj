// Action_Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Outtraj::Action_Outtraj() :
  max_(0.0),
  min_(0.0),
  Dset_(NULL)
{
  //fprintf(stderr,"Outtraj Con\n");
} 

// Action_Outtraj::init()
/** Expected call: outtraj <filename> [ trajout args ] 
  *                        [maxmin <dataset> min <min> max <max>
  */
int Action_Outtraj::init() {
#ifdef MPI
  mprintf("ERROR: OUTTRAJ currently not functional with MPI.\n");
  return 1;
#endif

  outtraj_.SetDebug(debug);
  Topology* tempParm = PFL->GetParm(actionArgs);
  if (tempParm==NULL) {
    mprinterr("Error: OUTTRAJ: Could not get parm for %s\n",actionArgs.ArgAt(1));
    return 1;
  }
  if ( outtraj_.SetupTrajWrite(actionArgs.GetStringNext(), &actionArgs, 
                               tempParm, TrajectoryFile::UNKNOWN_TRAJ) ) 
    return 1;
  mprintf("    OUTTRAJ:");
  outtraj_.PrintInfo(1);
  // If maxmin, get the name of the dataset as well as the max and min values.
  ArgList::ConstArg datasetName = actionArgs.getKeyString("maxmin");
  if (datasetName!=NULL) {
    Dset_ = DSL->Get(datasetName);
    if (Dset_==NULL) {
      mprintf("Error: Outtraj maxmin: Could not get dataset %s\n",datasetName);
      return 1;
    } else {
      // Currently only allow int, float, or double datasets
      if (Dset_->Type() != DataSet::INT &&
          Dset_->Type() != DataSet::FLOAT &&
          Dset_->Type() != DataSet::DOUBLE) 
      {
        mprinterr("Error: Outtraj maxmin: Only int, float, or double dataset (%s) supported.\n",
                datasetName);
        return 1;
      }
      max_ = actionArgs.getKeyDouble("max",0.0);
      min_ = actionArgs.getKeyDouble("min",0.0);
      mprintf("             maxmin: Printing trajectory frames based on %lf <= %s <= %lf\n",
              min_, datasetName, max_);
    }
  }

  return 0;
} 

// Action_Outtraj::action()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
int Action_Outtraj::action() {
  // If dataset defined, check if frame is within max/min
  if (Dset_!=NULL) {
    double dVal = Dset_->CurrentDval();
    //mprintf("DBG: maxmin: dVal = %lf\n",dVal);
    // If value from dataset not within min/max, exit now.
    if (dVal < min_ || dVal > max_) return 0;
  }
  if ( outtraj_.WriteFrame(frameNum, currentParm, *currentFrame) != 0 ) 
    return 1;
  return 0;
} 

// Action_Outtraj::print()
/** Close trajectory. Indicate how many frames were actually written.
  */
void Action_Outtraj::print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",actionArgs.ArgAt(1),outtraj_.NumFramesProcessed());
  outtraj_.EndTraj();
}

