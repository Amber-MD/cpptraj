#include "EnsembleOut_Multi.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // AppendNumber
#ifdef MPI
# include "Parallel.h"
#endif

// CONSTRUCTOR
EnsembleOut_Multi::EnsembleOut_Multi() : ensembleSize_(0) {} 

// DESTRUCTOR
EnsembleOut_Multi::~EnsembleOut_Multi() {
  EndEnsemble();
  Clear();
}

void EnsembleOut_Multi::Clear() {
  for (IOarrayType::const_iterator tio = ioarray_.begin(); tio != ioarray_.end(); ++tio)
    delete *tio;
  ioarray_.clear();
}

// EnsembleOut_Multi::InitEnsembleWrite()
/** Initialize each output trajectory with appropriate TrajectoryIO class
  * and process arguments.
  */
int EnsembleOut_Multi::InitEnsembleWrite(std::string const& tnameIn,
                                         ArgList const& argIn, int ensembleSizeIn,
                                         TrajectoryFile::TrajFormatType writeFormatIn)
{
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  ensembleSize_ = ensembleSizeIn; 
  // Require that ensemble size is set.
  if (ensembleSize_ < 1) {
    mprinterr("Internal Error: Ensemble size too small for ensemble output.\n");
    return 1;
  }
  ArgList trajout_args = argIn;
  // Get onlymembers range
  if (SetMembersToWrite(trajout_args.GetStringKey("onlymembers"), ensembleSize_)) return 1;
  // Process common args
  if (SetTraj().CommonTrajoutSetup(tnameIn, trajout_args, writeFormatIn))
    return 1;
  Clear();
  // Set up ensemble file names.
  fileNames_.clear();
# ifdef MPI
  // In MPI each thread writes a single member.
  if (MembersToWrite().InRange( Parallel::EnsembleComm().Rank() ))
    fileNames_.push_back( AppendNumber(Traj().Filename().Full(), Parallel::EnsembleComm().Rank()) );
  //else
  //  rprintf("DEBUG: Skipping member '%s'\n", 
  //          AppendNumber(Traj().Filename().Full(), Parallel::EnsembleComm().Rank()).c_str());
# else
  // In serial single process writes each member.
  // Create a map: tIndex[ pos ] = <ioarray_index>
  tIndex_.clear();
  tIndex_.reserve( ensembleSize_ );
  int ioidx = 0;
  for (int num = 0; num < ensembleSize_; num++) {
    if (MembersToWrite().InRange( num )) {
      fileNames_.push_back( AppendNumber(Traj().Filename().Full(), num) );
      tIndex_.push_back( ioidx++ );
    } else {
      //mprintf("DEBUG: Skipping member '%s'\n",
      //        AppendNumber(Traj().Filename().Full(), num).c_str());
      tIndex_.push_back( -1 );
    }
  }
# endif
  // Set up write format for each file. 
  typedef std::vector<TrajectoryFile::TrajFormatType> FmtArray;
  FmtArray fileFormats(fileNames_.size(), Traj().WriteFormat());
  // If appending, all files must exist and must have same format.
  if (Traj().Append()) {
    for (unsigned int m = 0; m != fileNames_.size(); ++m) {
      if (Traj().CheckAppendFormat( fileNames_[m], fileFormats[m] )) {
        mprintf("Warning: 'append' disabled; must be valid for all ensemble members.\n");
        // TODO: OnlyMembers-aware?
        SetTraj().SetAppend( false );
        break;
      }
    }
  }
  // Set up TrajectoryIO for each member.
  TrajectoryFile::TrajFormatType lastFormat = TrajectoryFile::UNKNOWN_TRAJ;
  for (unsigned int m = 0; m != fileNames_.size(); ++m) {
    if (fileFormats[m] != lastFormat)
      mprintf("\tWriting ensemble member '%s' as %s\n", fileNames_[m].c_str(),
              TrajectoryFile::FormatString(fileFormats[m]));
    lastFormat = fileFormats[m];
    TrajectoryIO* tio = TrajectoryFile::AllocTrajIO( fileFormats[m] );
    if (tio == 0) return 1;
    ioarray_.push_back( tio );
    ioarray_.back()->SetDebug( debug_ );
    // Set specified title - will not set if empty 
    //if (!TrajoutTitle().empty())
    //  trajio_->SetTitle( TrajoutTitle() + "." + integerToString( num ) );
    ioarray_.back()->SetTitle( Traj().Title() );
    // Process any write arguments specific to certain formats not related
    // to parm file. Options related to parm file are handled in SetupTrajWrite 
    ArgList rep_args = trajout_args;
    if (ioarray_.back()->processWriteArgs( rep_args )) {
      mprinterr("Error: trajout %s: Could not process arguments.\n",fileNames_[m].c_str());
      return 1;
    }
  }
  // Write is set up for topology in SetupTrajWrite 
  return 0;
}

// EnsembleOut_Multi::EndEnsemble()
void EnsembleOut_Multi::EndEnsemble() {
  //if (TrajIsOpen()) {
    for (IOarrayType::const_iterator tio = ioarray_.begin(); tio != ioarray_.end(); ++tio)
#     ifdef MPI
      if (trajComm_.Size() > 1)
        (*tio)->parallelCloseTraj();
      else
#     endif
        (*tio)->closeTraj();
  //  SetTrajIsOpen( false );
  //}
}

/** Perform any topology-related setup for this trajectory if given Topology
  * matches what trajectory was initialized with; the topology may have
  * been modified (e.g. by a 'strip' command) since the output trajectory was
  * initialized.
  */
int EnsembleOut_Multi::SetupEnsembleWrite(Topology* tparmIn, CoordinateInfo const& cInfoIn,
                                          int nFrames)
{
  // Setup topology and coordinate info.
  if (SetTraj().SetupCoordInfo(tparmIn, nFrames, cInfoIn))
    return 1;
# ifdef MPI
  if (!trajComm_.IsNull() && trajComm_.Size() > 1)
    return ParallelSetupEnsembleWrite();
# endif
  // Set up all TrajectoryIOs
  //if (!TrajIsOpen()) {
    for (unsigned int m = 0; m != ioarray_.size(); ++m) {
      if (ioarray_[m]->setupTrajout(fileNames_[m], Traj().Parm(), Traj().CoordInfo(),
                                    Traj().NframesToWrite(), Traj().Append()))
       return 1;
    }
  //}
  if (debug_ > 0)
    Traj().CoordInfo().PrintCoordInfo(Traj().Filename().base(), Traj().Parm()->c_str());
  return 0;
}

// EnsembleOut_Multi::WriteEnsemble()
/** Write given array of frames if trajectory is open (initialzed and set-up).
  */ 
int EnsembleOut_Multi::WriteEnsemble(int set, FramePtrArray const& Farray)
{
  // Check that set should be written
  if (SetTraj().CheckFrameRange(set)) return 0;
  // Write
# ifdef MPI
  int err = 0;
  if (!ioarray_.empty()) {
    if (trajComm_.Size() > 1)
      err = ioarray_.front()->parallelWriteFrame(set, *Farray.front());
    else
      err = ioarray_.front()->writeFrame(set, *Farray.front());
  }
  return err;
# else
  for (int member = 0; member != ensembleSize_; member++) {
    int tidx = tIndex_[member];
    if (tidx != -1) {
      if (ioarray_[tidx]->writeFrame(set, *Farray[member])) return 1;
    }
  }
  return 0;
# endif
}

// EnsembleOut_Multi::PrintInfo()
void EnsembleOut_Multi::PrintInfo(int expectedNframes) const {
  mprintf("  '%s.X' ", Traj().Filename().base());
  if (expectedNframes > 0) mprintf("(%i frames) ", expectedNframes);
  // Use MembersToWrite() instead of ioarray since 'onlymembers' may have been specified
  mprintf("(Ensemble, %i members written", MembersToWrite().Size());
  if (MembersToWrite().Size() < ensembleSize_) {
    mprintf(":");
    for (Range::const_iterator it = MembersToWrite().begin();
                               it != MembersToWrite().end(); ++it)
      mprintf(" %i", *it);
  }
  mprintf(") ");
  if (!ioarray_.empty()) //FIXME: If 'onlymembers' does not include 0, no info printed
    ioarray_.front()->Info();
  Traj().CommonInfo(); // NOTE: Prints newline
}
#ifdef MPI
// -----------------------------------------------------------------------------
int EnsembleOut_Multi::ParallelSetupEnsembleWrite()
{
  int err = 0;
  // Set up all TrajectoryIOs
  //if (!TrajIsOpen()) {
    for (unsigned int m = 0; m != ioarray_.size(); ++m) {
      err += ioarray_[m]->parallelSetupTrajout(fileNames_[m], Traj().Parm(), Traj().CoordInfo(),
                                               Traj().NframesToWrite(), Traj().Append(),
                                               trajComm_);
      err += ioarray_[m]->parallelOpenTrajout( trajComm_ );
    }
    if (trajComm_.CheckError( err )) return 1;
  //}
  return 0;
}
#endif
