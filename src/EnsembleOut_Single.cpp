#ifdef ENABLE_SINGLE_ENSEMBLE
#include "EnsembleOut_Single.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
EnsembleOut_Single::EnsembleOut_Single() : eio_(0), ensembleSize_(0) {}

// DESTRUCTOR
EnsembleOut_Single::~EnsembleOut_Single() {
  EndEnsemble();
  if (eio_ != 0) delete eio_;
}

// EnsembleOut_Single::InitEnsembleWrite()
int EnsembleOut_Single::InitEnsembleWrite(std::string const& tnameIn,
                                          ArgList const& argIn, int ensembleSizeIn,
                                          TrajectoryFile::TrajFormatType writeFormatIn)
{
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  // Require that ensemble size is set.
  ensembleSize_ = ensembleSizeIn; 
  if (ensembleSize_ < 1) {
    mprinterr("Internal Error: Ensemble size has not been set.\n");
    return 1;
  }
  ArgList trajout_args = argIn;
  // Get onlymembers range
  Range members_to_write = MembersToWrite(trajout_args.GetStringKey("onlymembers"), ensembleSize_);
  if (members_to_write.Empty()) return 1;
  // Process common args
  if (SetTraj().CommonTrajoutSetup(tnameIn, trajout_args, writeFormatIn))
    return 1;
  if (eio_ != 0) delete eio_;
  // If appending, file must exist and must match the current format.
  // TODO Do not pass in writeformat directly to be changed.
  if (Traj().Append()) {
    if (Traj().CheckAppendFormat( Traj().Filename(), SetTraj().WriteFormat() ))
      SetTraj().SetAppend( false );
  }
  // Set up for the specified format.
  eio_ = TrajectoryFile::AllocTrajIO( Traj().WriteFormat() );
  if (eio_ == 0) return 1;
  // Check that the TrajectoryIO object can read/write single ensemble
  if (!eio_->CanProcessEnsemble()) {
    mprinterr("Error: Format '%s' cannot be used for ensemble single file output.\n",
              TrajectoryFile::FormatString(Traj().WriteFormat()));
    return 1;
  }
  mprintf("\tWriting '%s' as %s\n", Traj().Filename().full(),
          TrajectoryFile::FormatString(Traj().WriteFormat()));
  eio_->SetDebug( debug_ );
  // Set specified title - will not set if empty 
  eio_->SetTitle( Traj().Title() );
  // Process any write arguments specific to certain formats not related
  // to parm file. Options related to parm file are handled in SetupTrajWrite 
  if (eio_->processWriteArgs(trajout_args)) {
    mprinterr("Error: trajout %s: Could not process arguments.\n",Traj().Filename().full());
    return 1;
  }
  return 0;
}

// EnsembleOut_Single::EndEnsemble()
void EnsembleOut_Single::EndEnsemble() {
  //if (TrajIsOpen()) {
  if (eio_ != 0) eio_->closeTraj(); // Handle no init case
  //  SetTrajIsOpen(false);
  //}
}

// EnsembleOut_Single::SetupEnsembleWrite()
int EnsembleOut_Single::SetupEnsembleWrite(Topology* tparmIn, CoordinateInfo const& cInfoIn,
                                           int nFrames)
{
  // Set up topology and coordinate info.
  if (SetTraj().SetupCoordInfo(tparmIn, nFrames, cInfoIn))
    return 1;
# ifdef MPI
  if (!trajComm_.IsNull() && trajComm_.Size() > 1)
    return ParallelSetupEnsembleWrite();
# endif
  if (debug_ > 0)
    rprintf("\tSetting up single ensemble %s for WRITE, topology '%s' (%i atoms).\n",
            Traj().Filename().base(), tparmIn->c_str(), tparmIn->Natom());
  // Set up TrajectoryIO
  if (eio_->setupTrajout(Traj().Filename(), Traj().Parm(), Traj().CoordInfo(),
                         Traj().NframesToWrite(), Traj().Append()))
    return 1;
  if (debug_ > 0)
    eio_->CoordInfo().PrintCoordInfo(Traj().Filename().base(), Traj().Parm()->c_str());
  // First frame setup
  //if (!TrajIsOpen()) { //}
  return 0;
}

// EnsembleOut_Single::WriteEnsemble()
int EnsembleOut_Single::WriteEnsemble(int set, FramePtrArray const& Farray) {
  // Check that set should be written
  if (SetTraj().CheckFrameRange(set)) return 0;
  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (eio_->writeArray(set, Farray)) return 1;
  return 0;
}

// EnsembleOut_Single::PrintInfo()
void EnsembleOut_Single::PrintInfo(int expectedNframes) const {
  mprintf("  '%s' (Single Ensemble, %i members",Traj().Filename().base(), ensembleSize_);
  if (expectedNframes > 0) mprintf(", %i frames", expectedNframes);
  mprintf(") ");
  eio_->Info();
  Traj().CommonInfo();
}
#ifdef MPI
int EnsembleOut_Single::ParallelSetupEnsembleWrite() {
  mprinterr("Error: Multiple threads per ensemble not supported for single ensemble write.\n");
  return 1;
}
#endif
#endif /* ENABLE_SINGLE_ENSEMBLE */
