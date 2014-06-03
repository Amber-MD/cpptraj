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

// Trajout_Single::InitStdoutTrajWrite()
int Trajout_Single::InitStdoutTrajWrite(ArgList const& argIn, Topology *tparmIn,
                                 TrajFormatType writeFormatIn)
{
  return InitTrajout("", argIn, tparmIn, writeFormatIn);
}

// Trajout_Single::InitEnsembleTrajWrite()
int Trajout_Single::InitEnsembleTrajWrite(std::string const& tnameIn, ArgList const& argIn,
                                          Topology* tparmIn, TrajFormatType fmtIn,
                                          int ensembleNum)
{
  FileName tempName;
  tempName.SetFileName( tnameIn );
  TrajFormatType extFmt = TrajectoryFile::GetTypeFromExtension( tempName.Ext() );
  if (extFmt != UNKNOWN_TRAJ)
    fmtIn = extFmt;
  if (ensembleNum > -1)
    return InitTrajWrite( NumberFilename(tnameIn, ensembleNum), argIn, tparmIn, extFmt );
  else
    return InitTrajWrite( tnameIn, argIn, tparmIn, extFmt );
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
  // to parm file. Options related to parm file are handled on the first
  // write in WriteFrame.
  if (trajio_->processWriteArgs(trajout_args)) {
    mprinterr("Error: trajout %s: Could not process arguments.\n",TrajFilename().full());
    return 1;
  }
  // Write is set up for topology only when first frame written.
  return 0;
}

// Trajout_Single::EndTraj()
void Trajout_Single::EndTraj() {
  if (TrajIsOpen()) {
    trajio_->closeTraj();
    SetTrajIsOpen(false);
  }
}

// Trajout_Single::WriteFrame()
/** Write given frame. If this is the first frame being written this routine
  * is where the output trajectory will actually be set up for the associated 
  * topology file since the topology may have been modified (e.g. by a 'strip'
  * command) since the output trajectory was initialized (modified topologies
  * will still have the same Pindex).
  */ 
int Trajout_Single::WriteFrame(int set, Topology *tparmIn, Frame const& FrameOut) {
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->Pindex() != TrajParm()->Pindex()) return 0;

  // First frame setup
  if (!TrajIsOpen()) {
    if (FirstFrameSetup(TrajFilename().Full(), trajio_, tparmIn)) return 1;
  }

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
