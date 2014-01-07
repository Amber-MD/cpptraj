// Action_Outtraj 
#include "Action_Outtraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Outtraj::Action_Outtraj() : CurrentParm_(0) {} 

void Action_Outtraj::Help() {
  mprintf("\t<filename> [ trajout args ]\n"
          "  Like 'trajout', but coordinate output occurs during actions rather than at the end.\n");
}

// Action_Outtraj::Init()
Action::RetType Action_Outtraj::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
# ifdef MPI
  mprintf("ERROR: OUTTRAJ currently not functional with MPI.\n");
  return Action::ERR;
# else
  // maxmin now deprecated, functionality is in Action_FilterByData
  if (actionArgs.Contains("maxmin") || actionArgs.Contains("max") ||
      actionArgs.Contains("min") || actionArgs.Contains("maxmindata")) {
    mprinterr("Error: maxmin functionality is deprecated. Please use the 'filter' action instead.\n");
    return Action::ERR;
  }
  // Set up output traj
  outtraj_.SetDebug(debugIn);
  std::string trajfilename = actionArgs.GetStringNext();
  if (trajfilename.empty()) {
    mprinterr("Error: No filename given.\nError: Usage: ");
    Help();
    return Action::ERR;
  }
  Topology* tempParm = PFL->GetParm(actionArgs);
  if (tempParm==0) {
    mprinterr("Error: OUTTRAJ: Could not get parm for %s\n",trajfilename.c_str());
    return Action::ERR;
  }
  if ( outtraj_.InitTrajWrite(trajfilename, actionArgs, 
                               tempParm, TrajectoryFile::UNKNOWN_TRAJ) ) 
    return Action::ERR;
  mprintf("    OUTTRAJ:");
  outtraj_.PrintInfo(1);

  return Action::OK;
# endif
} 

// Action_Outtraj::Setup()
Action::RetType Action_Outtraj::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm; 
  return Action::OK;
}

// Action_Outtraj::DoAction()
/** If a dataset was specified for maxmin, check if this structure
  * satisfies the criteria; if so, write. Otherwise just write.
  */
Action::RetType Action_Outtraj::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  if ( outtraj_.WriteFrame(frameNum, CurrentParm_, *currentFrame) != 0 ) 
    return Action::ERR;
  return Action::OK;
}

// Action_Outtraj::Print()
/** Close trajectory. Indicate how many frames were actually written.
  */
void Action_Outtraj::Print() {
  mprintf("  OUTTRAJ: [%s] Wrote %i frames.\n",outtraj_.TrajFilename().base(),
          outtraj_.NumFramesProcessed());
  outtraj_.EndTraj();
}
