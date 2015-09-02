#include "Action_LESsplit.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename

// DESTRUCTOR
Action_LESsplit::~Action_LESsplit() { if (lesParm_ != 0) delete lesParm_; }

void Action_LESsplit::Help() {
  mprintf("\t[out <filename prefix>] [average <avg filename>] <trajout args>\n"
          "  Split and/or average LES trajectory. At least one of 'out' or 'average'\n"
          "  must be specified. If both are specified they share <trajout args>.\n");
}

// Action_LESsplit::Init()
Action::RetType Action_LESsplit::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  if (DSL->EnsembleNum() > -1) {
    mprinterr("Error: LESSPLIT currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  trajfilename_ = actionArgs.GetStringKey("out");
  avgfilename_ = actionArgs.GetStringKey("average");
  lesSplit_ = !trajfilename_.empty();
  lesAverage_ = !avgfilename_.empty();
  if (!lesSplit_ && !lesAverage_) {
    mprinterr("Error: Must specify at least 'out <prefix>' or 'average <name>'.\n");
    return Action::ERR;
  }
  trajArgs_ = actionArgs.RemainingArgs();
  
  mprintf("    LESSPLIT:\n");
  if (lesSplit_) mprintf("\tSplit output to '%s.X'\n", trajfilename_.c_str());
  if (lesAverage_) mprintf("\tAverage output to '%s'\n", avgfilename_.c_str());
  return Action::OK;
}

// Action_LESsplit::Setup()
Action::RetType Action_LESsplit::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( !currentParm->LES().HasLES() ) {
    mprintf("Warning: No LES parameters in '%s', skipping.\n", currentParm->c_str());
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
    // Set up frames to hold individual copies
    lesFrames_.resize( lesMasks_.size() );
    lesFrames_.SetupFrames(lesParm_->Atoms(), lesParm_->ParmCoordInfo());
    lesPtrs_.resize( lesMasks_.size() );
    for (unsigned int i = 0; i != lesMasks_.size(); i++)
      lesPtrs_[i] = &lesFrames_[i];
    if (lesSplit_) {
      // Set up output ensemble FIXME check overwrites TODO combine init/setup?
      if (lesTraj_.InitEnsembleWrite(trajfilename_, trajArgs_, lesMasks_.size(),
                                     TrajectoryFile::UNKNOWN_TRAJ))
        return Action::ERR;
      if (lesTraj_.SetupEnsembleWrite(lesParm_, lesParm_->ParmCoordInfo(), lesParm_->Nframes()))
         return Action::ERR;
      lesTraj_.PrintInfo( 1 );
    }
    if (lesAverage_) {
      // For average only care about coords.
      avgFrame_.SetupFrame( lesParm_->Natom() );
      if (avgTraj_.PrepareTrajWrite( avgfilename_, trajArgs_, lesParm_,
                                     lesParm_->ParmCoordInfo(), lesParm_->Nframes(),
                                     TrajectoryFile::UNKNOWN_TRAJ ))
        return Action::ERR;
      avgTraj_.PrintInfo(1);
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
  for (unsigned int i = 0; i != lesMasks_.size(); i++)
    lesFrames_[i].SetFrame(*currentFrame, lesMasks_[i]);
  if (lesSplit_) {
    if ( lesTraj_.WriteEnsemble(frameNum, lesPtrs_) ) return Action::ERR;
  }
  if (lesAverage_) {
    avgFrame_.ZeroCoords();
    for (unsigned int i = 0; i != lesMasks_.size(); i++)
      avgFrame_ += lesFrames_[i];
    avgFrame_.Divide( lesMasks_.size() );
    if ( avgTraj_.WriteSingle(frameNum, avgFrame_) != 0 )
      return Action::ERR;
  }
  return Action::OK;
}
