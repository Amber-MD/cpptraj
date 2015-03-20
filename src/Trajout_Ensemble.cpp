#ifdef ENABLE_SINGLE_ENSEMBLE
#include "Trajout_Ensemble.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Trajout_Ensemble::Trajout_Ensemble() : eio_(0), ensembleSize_(0) {}

// DESTRUCTOR
Trajout_Ensemble::~Trajout_Ensemble() {
  EndTraj();
  if (eio_ != 0) delete eio_;
}

// Trajout_Ensemble::InitTrajWrite()
int Trajout_Ensemble::InitTrajWrite(std::string const& tnameIn, ArgList const& argIn,
                                 Topology *tparmIn, TrajFormatType writeFormatIn)
{
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  // Require that ensemble size is set.
  ensembleSize_ = tparmIn->ParmCoordInfo().EnsembleSize();
  if (ensembleSize_ < 1) {
    mprinterr("Internal Error: Ensemble size has not been set.\n");
    return 1;
  }
  ArgList trajout_args = argIn;
  TrajectoryFile::TrajFormatType writeFormat = writeFormatIn;
  // Get onlymembers range
  Range members_to_write = MembersToWrite(trajout_args.GetStringKey("onlymembers"), ensembleSize_);
  if (members_to_write.Empty()) return 1;
  // Process common args
  if (CommonTrajoutSetup(tnameIn, trajout_args, tparmIn, writeFormat))
    return 1;
  if (eio_ != 0) delete eio_;
  // If appending, file must exist and must match the current format.
  if (TrajoutAppend())
    CheckAppendFormat( tnameIn, writeFormat );
  // Set up for the specified format.
  eio_ = AllocTrajIO( writeFormat );
  if (eio_ == 0) return 1;
  // Check that the TrajectoryIO object can read/write single ensemble
  if (!eio_->CanProcessEnsemble()) {
    mprinterr("Error: Format '%s' cannot be used for ensemble single file output.\n",
              FormatString(writeFormat));
    return 1;
  }
  mprintf("\tWriting '%s' as %s\n", TrajFilename().full(),
          TrajectoryFile::FormatString(writeFormat));
  eio_->SetDebug( debug_ );
  // Set specified title - will not set if empty 
  eio_->SetTitle( TrajoutTitle() );
  // Process any write arguments specific to certain formats not related
  // to parm file. Options related to parm file are handled in SetupTrajWrite 
  if (eio_->processWriteArgs(trajout_args)) {
    mprinterr("Error: trajout %s: Could not process arguments.\n",TrajFilename().full());
    return 1;
  }
  return 0;
}

// Trajout_Ensemble::EndTraj()
void Trajout_Ensemble::EndTraj() {
  if (TrajIsOpen()) {
    eio_->closeTraj();
    SetTrajIsOpen(false);
  }
}

// Trajout_Ensemble::SetupTrajWrite()
int Trajout_Ensemble::SetupTrajWrite(Topology* tparmIn) {
  // First frame setup
  if (!TrajIsOpen()) {
    if (FirstFrameSetup(TrajFilename().Full(), eio_, tparmIn)) return 1;
  }
  return 0;
}

// Trajout_Ensemble::WriteEnsemble()
int Trajout_Ensemble::WriteEnsemble(int set, FramePtrArray const& Farray) {
  // Check that set should be written
  if (CheckFrameRange(set)) return 0;
  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (eio_->writeArray(set, Farray)) return 1;
  return 0;
}

// Trajout_Ensemble::PrintInfo()
void Trajout_Ensemble::PrintInfo(int showExtended) const {
  mprintf("  '%s' (Single Ensemble, %i members) ",TrajFilename().base(), ensembleSize_);
  CommonInfo( eio_ );
}
#endif
