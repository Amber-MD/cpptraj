#include "Action_LESsplit.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// DESTRUCTOR
Action_LESsplit::~Action_LESsplit() {
  for (TrajoutArray::iterator tout = lesTraj_.begin(); tout != lesTraj_.end(); ++tout)
    delete *tout;
}

void Action_LESsplit::Help() {
  mprintf("\tout <filename prefix> <trajout args>\n");
}

// Action_LESsplit::Init()
Action::RetType Action_LESsplit::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  trajfilename_ = actionArgs.GetStringKey("out");
  if (trajfilename_.empty()) {
    mprinterr("Error: Must specify output traj name prefix ('out <prefix>')\n");
    return Action::ERR;
  }
  trajArgs_ = actionArgs.RemainingArgs();
  
  mprintf("    LESSPLIT: Output to '%s.X'\n", trajfilename_.c_str());
  return Action::OK;
}

// Action_LESsplit::Setup()
Action::RetType Action_LESsplit::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->LES().Ntypes() < 1 ) {
    mprintf("Warning: No LES parameters in '%s', skipping.\n");
    return Action::ERR;
  }
  if (lesParm_ == 0) { // First time setup
    // Set up masks for all copies
    lesMasks_.clear();
    lesMasks_.resize( currentParm->LES().Ncopies() );
    unsigned int atom = 0;
    for (LES_Array::const_iterator les = currentParm->LES().Array().begin();
                                   les != currentParm->LES().Array().end(); ++les, ++atom)
    {
      // Copy 0 is in all copies
      if ( les->Copy() == 0 ) {
        for (MaskArray::iterator mask = lesMasks_.begin(); mask != lesMasks_.end(); ++mask)
          mask->AddAtom( atom );
      } else
        lesMasks_[ les->Copy() - 1 ].AddAtom( atom );
    }
    for (unsigned int i = 0; i < lesMasks_.size(); i++) {
      mprintf("\t%i atoms in LES copy %u\n", lesMasks_[i].Nselected(), i+1);
      if ( lesMasks_[i].Nselected() != lesMasks_[0].Nselected() ) {
        mprinterr("Error: Currently all LES copies MUST have same # atoms.\n");
        return Action::ERR;
      }
    }
    // Create topology for first copy
    lesParm_ = currentParm->modifyStateByMask( lesMasks_[0] );
    if (lesParm_ == 0) return Action::ERR;
    // Set up frame to hold individual copy
    lesFrame_.SetupFrameV(lesParm_->Atoms(), lesParm_->HasVelInfo(), lesParm_->NrepDim());
    // Set up output trajectories
    lesTraj_.clear();
    lesTraj_.reserve( lesMasks_.size() );
    for (unsigned int i = 0; i < lesMasks_.size(); i++) {
      lesTraj_.push_back( new Trajout() );
      // Copy trajArgs so they are the same for each.
      // FIXME: Should InitTrajWrite take const?
      ArgList targ = trajArgs_;
      if ( lesTraj_.back()->InitTrajWrite(NumberFilename( trajfilename_, i+1 ), targ,
                                          lesParm_, TrajectoryFile::UNKNOWN_TRAJ) )
        return Action::ERR;
      lesTraj_.back()->PrintInfo(1);
    }
  } else {
    if (lesParm_->Pindex() != currentParm->Pindex()) {
      mprintf("Warning: Already set up for LES parm '%s'. Skipping '%s'\n",
              lesParm_->c_str(), currentParm->c_str());
      return Action::ERR;
    }
  }

  return Action::OK;
}

// Action_LESsplit::DoAction()
Action::RetType Action_LESsplit::DoAction(int frameNum, Frame* currentFrame, 
                                          Frame** frameAddress)
{
  for (unsigned int i = 0; i < lesMasks_.size(); i++) {
    lesFrame_.SetFrame(*currentFrame, lesMasks_[i]);
    if ( lesTraj_[i]->WriteFrame(frameNum, lesParm_, lesFrame_) != 0 )
      return Action::ERR;
  }
  return Action::OK;
}
