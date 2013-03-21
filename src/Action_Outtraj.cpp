// Action_Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Outtraj::Action_Outtraj() :
  maxmin_(0),
  CurrentParm_(0)
{ } 

void Action_Outtraj::Help() {
  mprintf("\t<filename> [ trajout args ]\n");
  mprintf("\t[maxmin <dataset> min <min> max <max>] ...\n");
  mprintf("\tmaxmindata <file>\n");
  mprintf("\tLike 'trajout', but coordinate output occurs during actions rather than at the end.\n");
}

// Action_Outtraj::init()
Action::RetType Action_Outtraj::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
#ifdef MPI
  mprintf("ERROR: OUTTRAJ currently not functional with MPI.\n");
  return Action::ERR;
#endif
  outtraj_.SetDebug(debugIn);
  std::string trajfilename = actionArgs.GetStringNext();
  if (trajfilename.empty()) {
    mprinterr("Error: outtraj: no filename given.\nError: Usage: ");
    Help();
    return Action::ERR;
  }
  Topology* tempParm = PFL->GetParm(actionArgs);
  if (tempParm==0) {
    mprinterr("Error: OUTTRAJ: Could not get parm for %s\n",trajfilename.c_str());
    return Action::ERR;
  }
  if ( outtraj_.SetupTrajWrite(trajfilename, &actionArgs, 
                               tempParm, TrajectoryFile::UNKNOWN_TRAJ) ) 
    return Action::ERR;
  mprintf("    OUTTRAJ:");
  outtraj_.PrintInfo(1);
  // If maxmin, get the name of the dataset as well as the max and min values.
  while ( actionArgs.Contains("maxmin") ) {
    std::string datasetName = actionArgs.GetStringKey("maxmin");
    if (!datasetName.empty()) {
      DataSet* dset = DSL->GetDataSet(datasetName);
      if (dset==0) {
        mprintf("Error: Outtraj maxmin: Could not get dataset %s\n",datasetName.c_str());
        return Action::ERR;
      } else {
        // Currently only allow int, float, or double datasets
        if (dset->Type() != DataSet::INT &&
            dset->Type() != DataSet::FLOAT &&
            dset->Type() != DataSet::DOUBLE) 
        {
          mprinterr("Error: Outtraj maxmin: Only int, float, or double dataset (%s) supported.\n",
                  datasetName.c_str());
          return Action::ERR;
        }
        Dsets_.push_back( dset );
        Max_.push_back( actionArgs.getKeyDouble("max",0.0) );
        Min_.push_back( actionArgs.getKeyDouble("min",0.0) );
        mprintf("             maxmin: Printing trajectory frames based on %f <= %s <= %f\n",
                Min_.back(), datasetName.c_str(), Max_.back());
      }
    } else {
      mprinterr("Error: outtraj: maxmin Usage: maxmin <setname> <max> <min>\n");
      return Action::ERR;
    }
  }
  if (!Dsets_.empty()) {
    std::string maxmindata = actionArgs.GetStringKey("maxmindata");
    if (!maxmindata.empty()) {
      maxmin_ = DSL->AddSet( DataSet::INT, actionArgs.GetStringNext(), "maxmin" );
      if (maxmin_ == 0) return Action::ERR;
      DFL->AddSetToFile( maxmindata, maxmin_ );
      mprintf("\tMaxMin frame info will be written to %s\n", maxmindata.c_str());
    }
  }
  return Action::OK;
} 

Action::RetType Action_Outtraj::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm; 
  return Action::OK;
}

// Action_Outtraj::action()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
Action::RetType Action_Outtraj::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  static int ONE = 1;
  // If dataset defined, check if frame is within max/min
  if (!Dsets_.empty()) {
    for (unsigned int ds = 0; ds < Dsets_.size(); ++ds)
    {
      double dVal = Dsets_[ds]->Dval(frameNum);
      //mprintf("DBG: maxmin[%u]: dVal = %f, min = %f, max = %f\n",ds,dVal,Min_[ds],Max_[ds]);
      // If value from dataset not within min/max, exit now.
      if (dVal < Min_[ds] || dVal > Max_[ds]) return Action::OK;
    }
    if (maxmin_ != 0)
      maxmin_->Add( frameNum, &ONE );
  }
  if ( outtraj_.WriteFrame(frameNum, CurrentParm_, *currentFrame) != 0 ) 
    return Action::ERR;
  return Action::OK;
} 

// Action_Outtraj::print()
/** Close trajectory. Indicate how many frames were actually written.
  */
void Action_Outtraj::Print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",outtraj_.TrajFilename().base(),
          outtraj_.NumFramesProcessed());
  outtraj_.EndTraj();
}

