#include "Trajout_Single.h"
#include "TrajectoryIO.h"
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // AppendNumber
#include "ArgList.h"

// DESTRUCTOR
Trajout_Single::~Trajout_Single() {
  EndTraj();
  if (trajio_!=0) delete trajio_;
}

// Trajout_Single::InitTrajWrite()
/** Initialize output trajectory with appropriate TrajectoryIO class and 
  * process arguments.
  */
int Trajout_Single::InitTrajWrite(FileName const& tnameIn, ArgList const& argIn,
                                  DataSetList const& DSLin,
                                  TrajectoryFile::TrajFormatType writeFormatIn)
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  return InitTrajout(tnameIn, argIn, DSLin, writeFormatIn);
}

/** Initialize trajectory for write. Append ensemble number to filename if given. */
int Trajout_Single::InitEnsembleTrajWrite(FileName const& tnameIn, ArgList const& argIn,
                                          DataSetList const& DSLin, 
                                          TrajectoryFile::TrajFormatType fmtIn, int ensembleNum)
{
  ArgList args = argIn;
  // Try to determine format here, before extension is potentially modified.
  TrajectoryFile::TrajFormatType fmt = fmtIn;
  if (fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromArg( args, TrajectoryFile::UNKNOWN_TRAJ );
  if (fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromFname( tnameIn, TrajectoryFile::UNKNOWN_TRAJ );
  int err = 0;
  if (ensembleNum > -1)
    err = InitTrajWrite( AppendNumber(tnameIn.Full(), ensembleNum), args, DSLin, fmt );
  else
    err = InitTrajWrite( tnameIn, args, DSLin, fmt );
  if (err != 0) return 1;
  return 0;
}

// Trajout_Single::PrepareTrajWrite()
int Trajout_Single::PrepareTrajWrite(FileName const& tnameIn, ArgList const& argIn,
                                     DataSetList const& DSLin,
                                     Topology* tparmIn, CoordinateInfo const& cInfoIn,
                                     int nFrames, TrajectoryFile::TrajFormatType fmtIn)
{
  if (InitTrajWrite(tnameIn, argIn, DSLin, fmtIn)) return 1;
  if (SetupTrajWrite(tparmIn, cInfoIn, nFrames)) return 1;
  return 0;
}

// Trajout_Single::InitTrajout()
int Trajout_Single::InitTrajout(FileName const& tnameIn, ArgList const& argIn,
                                DataSetList const& DSLin,
                                TrajectoryFile::TrajFormatType writeFormatIn)
{
  ArgList trajout_args = argIn;
  // Process common args
  if (traj_.CommonTrajoutSetup(tnameIn, trajout_args, writeFormatIn))
    return 1;
  if (trajio_ != 0) delete trajio_;
  // If appending, file must exist and must match the current format.
  // TODO Do not pass in writeformat directly to be changed.
  if (traj_.Append()) {
    if (traj_.CheckAppendFormat( traj_.Filename(), traj_.WriteFormat() ))
      traj_.SetAppend( false );
  }
  // Set up for the specified format.
  trajio_ = TrajectoryFile::AllocTrajIO( traj_.WriteFormat() );
  if (trajio_ == 0) return 1;
  mprintf("\tWriting '%s' as %s\n", traj_.Filename().full(),
          TrajectoryFile::FormatString(traj_.WriteFormat()));
  trajio_->SetDebug( debug_ );
  // Set specified title - will not set if empty 
  trajio_->SetTitle( traj_.Title() );
  // Process any write arguments specific to certain formats not related
  // to parm file. Options related to parm file are handled in SetupTrajWrite 
  if (trajio_->processWriteArgs(trajout_args, DSLin)) {
    mprinterr("Error: trajout %s: Could not process arguments.\n", traj_.Filename().full());
    return 1;
  }
  // Write is set up for topology in SetupTrajWrite 
  return 0;
}

// Trajout_Single::EndTraj()
void Trajout_Single::EndTraj() {
  //if (TrajIsOpen()) { // FIXME: Necessary?
  if (trajio_ != 0) {  // Handle no init case
#   ifdef MPI
    if (trajComm_.Size() > 1)
      trajio_->parallelCloseTraj();
    else
#   endif
      trajio_->closeTraj();
  }
  //  SetTrajIsOpen(false);
  //}
}

/** Perform any topology-related setup for this trajectory */
int Trajout_Single::SetupTrajWrite(Topology* tparmIn, CoordinateInfo const& cInfoIn, int nFrames) {
  // Set up topology and coordinate info.
  if (traj_.SetupCoordInfo(tparmIn, nFrames, cInfoIn))
    return 1;
# ifdef MPI
  if (!trajComm_.IsNull() && trajComm_.Size() > 1)
    return ParallelSetupTrajWrite();
# endif
  if (debug_ > 0)
    rprintf("\tSetting up %s for WRITE, topology '%s' (%i atoms).\n",
            traj_.Filename().base(), tparmIn->c_str(), tparmIn->Natom());
  // Set up TrajectoryIO
  if (trajio_->setupTrajout(traj_.Filename(), traj_.Parm(), traj_.CoordInfo(),
                            traj_.NframesToWrite(), traj_.Append()))
    return 1;
  if (debug_ > 0)
    trajio_->CoordInfo().PrintCoordInfo(traj_.Filename().base(), traj_.Parm()->c_str());
  // First frame setup
  //if (!TrajIsOpen()) { //}
  return 0;
}

// Trajout_Single::WriteSingle()
/** Write given frame if trajectory is open (initialzed and set-up).
  */
int Trajout_Single::WriteSingle(int set, Frame const& FrameOut) {
  // Check that set should be written
  if (traj_.CheckFrameRange(set)) return 0;
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  // Write
  int err = 0;
# ifdef MPI
  if (trajComm_.Size() > 1)
    err = trajio_->parallelWriteFrame(set, FrameOut);
  else
# endif
    err = trajio_->writeFrame(set, FrameOut);
  return err;
}

// Trajout_Single::PrintInfo()
void Trajout_Single::PrintInfo(int expectedNframes) const {
  mprintf("  '%s' ", traj_.Filename().base());
  if (expectedNframes > 0) mprintf("(%i frames) ", expectedNframes);
  trajio_->Info();
  traj_.CommonInfo();
}
#ifdef MPI
// -----------------------------------------------------------------------------
int Trajout_Single::ParallelSetupTrajWrite() {
  if (traj_.HasRange()) {
    mprinterr("Error: 'trajout' frame options not yet supported during parallel processing.\n");
    return 1;
  }
  if (traj_.Append()) {
    mprinterr("Error: 'trajout' append not yet supported during parallel processing.\n");
    return 1;
  }
  if (debug_ > 0)
    rprintf("\tSetting up '%s' for WRITE in parallel, topology '%s' (%i atoms, %i frames).\n",
            traj_.Filename().base(), traj_.Parm()->c_str(), traj_.Parm()->Natom(),
            traj_.NframesToWrite());
  // Set up TrajectoryIO in parallel.
  if (trajio_->parallelSetupTrajout(traj_.Filename(), traj_.Parm(), traj_.CoordInfo(),
                                    traj_.NframesToWrite(), traj_.Append(), trajComm_))
  {
    mprinterr("Error: Could not set up parallel trajout.\n");
    mprinterr("Error: Problem with set up or format may be unsupported in parallel.\n");
    return 1;
  }
  //if (debug_ > 0)
  //  Frame::PrintCoordInfo(traj_.Filename().base(), traj_.Parm()->c_str(), trajio_->CoordInfo());
  // Open TrajectoryIO in parallel.
  if (trajio_->parallelOpenTrajout( trajComm_ ))
  {
    mprinterr("Error: Could not open parallel trajout.\n");
    return 1;
  }
  return 0;
}
#endif
