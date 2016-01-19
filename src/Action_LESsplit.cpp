#include "Action_LESsplit.h"
#include "CpptrajStdio.h"

// DESTRUCTOR
Action_LESsplit::~Action_LESsplit() {
  if (lesSplit_) {
    for (Tarray::iterator tout = lesTraj_.begin(); tout != lesTraj_.end(); ++tout) {
      (*tout)->EndTraj();
      delete *tout;
    }
  }
  if (lesParm_ != 0) delete lesParm_;
}

void Action_LESsplit::Help() const {
  mprintf("\t[out <filename prefix>] [average <avg filename>] <trajout args>\n"
          "  Split and/or average LES trajectory. At least one of 'out' or 'average'\n"
          "  must be specified. If both are specified they share <trajout args>.\n");
}

// Action_LESsplit::Init()
Action::RetType Action_LESsplit::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  if (init.DSL().EnsembleNum() > -1) {
    mprinterr("Error: LESSPLIT currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  splitfilename_ = actionArgs.GetStringKey("out");
  std::string avgfilename = actionArgs.GetStringKey("average");
  lesSplit_ = !splitfilename_.empty();
  lesAverage_ = !avgfilename.empty();
  if (!lesSplit_ && !lesAverage_) {
    mprinterr("Error: Must specify at least 'out <prefix>' or 'average <name>'.\n");
    return Action::ERR;
  }
  trajArgs_ = actionArgs.RemainingArgs();
  // Initialize average output trajectory.
  // NOTE: Cannot yet init split traj since we dont know how many.
  if (lesAverage_) {
    avgTraj_.SetDebug( debugIn );
    if (avgTraj_.InitTrajWrite( avgfilename, trajArgs_, TrajectoryFile::UNKNOWN_TRAJ ))
      return Action::ERR;
  }
  
  mprintf("    LESSPLIT:\n");
  if (lesSplit_) mprintf("\tSplit output to '%s.X'\n", splitfilename_.c_str());
  if (lesAverage_) mprintf("\tAverage output to '%s'\n", avgTraj_.Traj().Filename().full());
  return Action::OK;
}

#ifdef MPI
int Action_LESsplit::ParallelActionInit(Parallel::Comm const& commIn) {
  trajComm_ = commIn;
  if (lesAverage_) avgTraj_.SetTrajComm( commIn );
  return 0;
}
#endif

// Action_LESsplit::Setup()
Action::RetType Action_LESsplit::Setup(ActionSetup& setup) {
  if ( !setup.Top().LES().HasLES() ) {
    mprintf("Warning: No LES parameters in '%s', skipping.\n", setup.Top().c_str());
    return Action::SKIP;
  }
  if (lesParm_ == 0) { // First time setup
    // Set up masks for all copies
    lesMasks_.clear();
    lesMasks_.resize( setup.Top().LES().Ncopies() );
    unsigned int atom = 0;
    for (LES_Array::const_iterator les = setup.Top().LES().Array().begin();
                                   les != setup.Top().LES().Array().end(); ++les, ++atom)
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
    lesParm_ = setup.Top().modifyStateByMask( lesMasks_[0] );
    if (lesParm_ == 0) return Action::ERR;
    // Set up frame to hold individual copy
    lesFrame_.SetupFrameV( lesParm_->Atoms(), setup.CoordInfo() );
    if (lesSplit_) {
      // Set up split output trajectories. lesTraj_ should be empty.
      lesTraj_.reserve( lesMasks_.size() );
      for (unsigned int idx = 0; idx != lesMasks_.size(); idx++) {
        // FIXME this will have to be changed if lessplit every enabled for ensemble
        lesTraj_.push_back( new Trajout_Single() );
        if (lesTraj_.back()->InitEnsembleTrajWrite( splitfilename_, trajArgs_,
                                                    TrajectoryFile::UNKNOWN_TRAJ, idx ))
          return Action::ERR;
#       ifdef MPI
        lesTraj_.back()->SetTrajComm( trajComm_ );
#       endif
        if (lesTraj_.back()->SetupTrajWrite( lesParm_, setup.CoordInfo(), setup.Nframes() ))
          return Action::ERR;
        lesTraj_.back()->PrintInfo(0);
      }
    }
    if (lesAverage_) {
      // For average only care about coords.
      avgFrame_.SetupFrame( lesParm_->Natom() );
      if (avgTraj_.SetupTrajWrite( lesParm_, CoordinateInfo(), setup.Nframes() ))
        return Action::ERR;
      avgTraj_.PrintInfo(0);
    }
  } else {
    if (lesParm_->Pindex() != setup.Top().Pindex()) {
      mprintf("Warning: Already set up for LES parm '%s'. Skipping '%s'\n",
              lesParm_->c_str(), setup.Top().c_str());
      return Action::SKIP;
    }
  }

  return Action::OK;
}

// Action_LESsplit::DoAction()
Action::RetType Action_LESsplit::DoAction(int frameNum, ActionFrame& frm) {
  if (lesAverage_)
    avgFrame_.ZeroCoords();
  for (unsigned int idx = 0; idx != lesMasks_.size(); idx++) {
    lesFrame_.SetFrame(frm.Frm(), lesMasks_[idx]);
    if (lesAverage_)
      avgFrame_ += lesFrame_;
    if (lesSplit_)
      if ( lesTraj_[idx]->WriteSingle(frm.TrajoutNum(), lesFrame_) != 0 )
        return Action::ERR;
  }
  if (lesAverage_) {
    avgFrame_.Divide( (double)lesMasks_.size() );
    if ( avgTraj_.WriteSingle(frm.TrajoutNum(), avgFrame_) != 0 )
      return Action::ERR;
  }
  return Action::OK;
}
