// Action_Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Outtraj::Action_Outtraj() :
  maxmin_(0)
{ } 

void Action_Outtraj::Help() {

}

// Action_Outtraj::init()
/** Expected call: outtraj <filename> [ trajout args ]
  *                        [maxmin <dataset> min <min> max <max>] ...
  *                        maxmindata <file>
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
  while ( actionArgs.Contains("maxmin") ) {
    ArgList::ConstArg datasetName = actionArgs.getKeyString("maxmin");
    if (datasetName!=NULL) {
      DataSet* dset = DSL->Get(datasetName);
      if (dset==NULL) {
        mprintf("Error: Outtraj maxmin: Could not get dataset %s\n",datasetName);
        return 1;
      } else {
        // Currently only allow int, float, or double datasets
        if (dset->Type() != DataSet::INT &&
            dset->Type() != DataSet::FLOAT &&
            dset->Type() != DataSet::DOUBLE) 
        {
          mprinterr("Error: Outtraj maxmin: Only int, float, or double dataset (%s) supported.\n",
                  datasetName);
          return 1;
        }
        Dsets_.push_back( dset );
        Max_.push_back( actionArgs.getKeyDouble("max",0.0) );
        Min_.push_back( actionArgs.getKeyDouble("min",0.0) );
        mprintf("             maxmin: Printing trajectory frames based on %f <= %s <= %f\n",
                Min_.back(), datasetName, Max_.back());
      }
    } else {
      mprinterr("Error: outtraj: maxmin Usage: maxmin <setname> <max> <min>\n");
      return 1;
    }
  }
  if (!Dsets_.empty()) {
    std::string maxmindata = actionArgs.GetStringKey("maxmindata");
    if (!maxmindata.empty()) {
      maxmin_ = DSL->AddSet( DataSet::INT, actionArgs.GetStringNext(), "maxmin" );
      if (maxmin_ == 0) return 1;
      DFL->AddSetToFile( maxmindata, maxmin_ );
      mprintf("\tMaxMin frame info will be written to %s\n", maxmindata.c_str());
    }
  }
  
  return 0;
} 

// Action_Outtraj::action()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
int Action_Outtraj::action() {
  static int ONE = 1;
  // If dataset defined, check if frame is within max/min
  if (!Dsets_.empty()) {
    for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    {
      double dVal = Dsets_[ds]->Dval(frameNum);
      //mprintf("DBG: maxmin[%u]: dVal = %f, min = %f, max = %f\n",ds,dVal,Min_[ds],Max_[ds]);
      // If value from dataset not within min/max, exit now.
      if (dVal < Min_[ds] || dVal > Max_[ds]) return 0;
    }
    if (maxmin_ != 0)
      maxmin_->Add( frameNum, &ONE );
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

