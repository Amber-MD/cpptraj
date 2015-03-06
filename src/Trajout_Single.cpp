#include "Trajout_Single.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" //NumberFilename

// CONSTRUCTOR
Trajout_Single::Trajout_Single() : trajio_(0) {}

// DESTRUCTOR
Trajout_Single::~Trajout_Single() {
  EndTraj();
  if (trajio_!=0) delete trajio_;
}

// Trajout_Single::InitTrajWrite()
/** Initialize output trajectory with appropriate TrajectoryIO class and 
  * process arguments.
  */
int Trajout_Single::InitTrajWrite(std::string const& tnameIn, ArgList const& argIn, 
                           Topology *tparmIn, TrajFormatType writeFormatIn)
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  return InitTrajout(tnameIn, argIn, tparmIn, writeFormatIn);
}

// Trajout_Single::PrepareStdoutTrajWrite()
/** Initialize and set up output trajectory for STDOUT write. */
int Trajout_Single::PrepareStdoutTrajWrite(ArgList const& argIn, Topology *tparmIn,
                                 TrajFormatType writeFormatIn)
{
  if (InitTrajout("", argIn, tparmIn, writeFormatIn)) return 1;
  if (SetupTrajWrite(tparmIn)) return 1;
  return 0;
}

// Trajout_Single::InitEnsembleTrajWrite()
int Trajout_Single::InitEnsembleTrajWrite(std::string const& tnameIn, ArgList const& argIn,
                                          Topology* tparmIn, TrajFormatType fmtIn,
                                          int ensembleNum)
{
  FileName tempName;
  tempName.SetFileName( tnameIn );
  TrajFormatType fmt = fmtIn;
  if (fmt == UNKNOWN_TRAJ)
    fmt = TrajectoryFile::GetTypeFromExtension( tempName.Ext() );
  if (ensembleNum > -1)
    return InitTrajWrite( NumberFilename(tnameIn, ensembleNum), argIn, tparmIn, fmt );
  else
    return InitTrajWrite( tnameIn, argIn, tparmIn, fmt );
}

// Trajout_Single::InitTrajout()
int Trajout_Single::InitTrajout(std::string const& tnameIn, ArgList const& argIn,
                                Topology *tparmIn, TrajFormatType writeFormatIn)
{
  ArgList trajout_args = argIn;
  TrajectoryFile::TrajFormatType writeFormat = writeFormatIn;
  // Process common args
  if (CommonTrajoutSetup(tnameIn, trajout_args, tparmIn, writeFormat))
    return 1;
  if (trajio_ != 0) delete trajio_;
  // If appending, file must exist and must match the current format.
  if (TrajoutAppend()) 
    CheckAppendFormat( tnameIn, writeFormat );
  // Set up for the specified format.
  trajio_ = AllocTrajIO( writeFormat );
  if (trajio_ == 0) return 1;
  mprintf("\tWriting '%s' as %s\n", TrajFilename().full(), 
          TrajectoryFile::FormatString(writeFormat));
  trajio_->SetDebug( debug_ );
  // Set specified title - will not set if empty 
  trajio_->SetTitle( TrajoutTitle() );
  // Process any write arguments specific to certain formats not related
  // to parm file. Options related to parm file are handled in SetupTrajWrite 
  if (trajio_->processWriteArgs(trajout_args)) {
    mprinterr("Error: trajout %s: Could not process arguments.\n",TrajFilename().full());
    return 1;
  }
  // Write is set up for topology in SetupTrajWrite 
  return 0;
}

// Trajout_Single::EndTraj()
void Trajout_Single::EndTraj() {
  if (TrajIsOpen()) {
    trajio_->closeTraj();
    SetTrajIsOpen(false);
  }
}

/** Perform any topology-related setup for this trajectory if given Topology
  * Pindex matches what trajectory was initialized with; the topology may have
  * been modified (e.g. by a 'strip' command) since the output trajectory was
  * initialized.
  */
int Trajout_Single::SetupTrajWrite(Topology* tparmIn) {
  // First frame setup
  if (!TrajIsOpen()) {
    if (FirstFrameSetup(TrajFilename().Full(), trajio_, tparmIn)) return 1;
  }
  return 0;
}

// Trajout_Single::WriteSingle()
/** Write given frame if trajectory is open (initialzed and set-up).
  */ 
int Trajout_Single::WriteSingle(int set, Frame const& FrameOut) {
  // Check that set should be written
  if (CheckFrameRange(set)) return 0;
  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (trajio_->writeFrame(set, FrameOut)) return 1;
  return 0;
}

// Trajout_Single::PrintInfo()
void Trajout_Single::PrintInfo(int showExtended) const {
  mprintf("  '%s' ",TrajFilename().base());
  CommonInfo( trajio_ );
}
