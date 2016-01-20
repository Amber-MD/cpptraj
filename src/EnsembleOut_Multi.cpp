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
  Range members_to_write = MembersToWrite(trajout_args.GetStringKey("onlymembers"), ensembleSize_);
  if (members_to_write.Empty()) return 1;
  // DEBUG
  //std::string dbg_mtw = "MembersToWrite:";
  //for (Range::const_iterator r = members_to_write.begin(); r != members_to_write.end(); r++)
  //  dbg_mtw += (" " + integerToString(*r));
  //rprintf("DEBUG: %s\n", dbg_mtw.c_str());
  // Process common args
  if (SetTraj().CommonTrajoutSetup(tnameIn, trajout_args, writeFormatIn))
    return 1;
  Clear();
  // Set up ensemble file names.
  fileNames_.clear();
# ifdef MPI
  // In MPI each thread writes a single member.
  if (members_to_write.InRange( Parallel::EnsembleComm().Rank() ))
    fileNames_.push_back( AppendNumber(Traj().Filename().Full(), Parallel::EnsembleComm().Rank()) );
  else
    rprintf("Warning: Skipping member '%s'\n", 
            AppendNumber(Traj().Filename().Full(), Parallel::EnsembleComm().Rank()).c_str());
# else
  // In serial single process writes each member.
  // Create a map: tIndex[ pos ] = <ioarray_index>
  tIndex_.clear();
  tIndex_.reserve( ensembleSize_ );
  int ioidx = 0;
  for (int num = 0; num < ensembleSize_; num++) {
    if (members_to_write.InRange( num )) {
      fileNames_.push_back( AppendNumber(Traj().Filename().Full(), num) );
      tIndex_.push_back( ioidx++ );
    } else {
      mprintf("Warning: Skipping member '%s'\n",
              AppendNumber(Traj().Filename().Full(), num).c_str());
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
  for (unsigned int m = 0; m != fileNames_.size(); ++m) {
    if (debug_ > 0)
      rprintf("\tWriting ensemble member '%s' as %s\n", fileNames_[m].c_str(),
              TrajectoryFile::FormatString(fileFormats[m]));
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
      (*tio)->closeTraj();
  //  SetTrajIsOpen( false );
  //}
}

/** Perform any topology-related setup for this trajectory if given Topology
  * matches what trajectory was initialized with; the topology may have
  * been modified (e.g. by a 'strip' command) since the output trajectory was
  * initialized.
  */
int EnsembleOut_Multi::SetupEnsembleWrite(Topology* tparmIn, CoordinateInfo const& cInfoIn, int nFrames) {
  // Setup topology and coordiante info.
  if (SetTraj().SetupCoordInfo(tparmIn, nFrames, cInfoIn))
    return 1;
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
  if (!ioarray_.empty()) {
    if (ioarray_.front()->writeFrame(set, *Farray.front())) return 1;
  }
# else
  for (int member = 0; member != ensembleSize_; member++) {
    int tidx = tIndex_[member];
    if (tidx != -1) {
      if (ioarray_[tidx]->writeFrame(set, *Farray[member])) return 1;
    }
  }
# endif
  return 0;
}

// EnsembleOut_Multi::PrintInfo()
void EnsembleOut_Multi::PrintInfo(int expectedNframes) const {
  mprintf("  '%s.X' ", Traj().Filename().base());
  if (expectedNframes > 0) mprintf("(%i frames) ", expectedNframes);
  mprintf("(Ensemble,");
# ifdef MPI
  // Since not every thread may be writing if 'onlymembers' specified,
  // determine total number being written.
  int mysize = (int)ioarray_.size();
  int total;
  Parallel::EnsembleComm().Reduce(&total, &mysize, 1, MPI_INT, MPI_SUM);
  mprintf(" %i members written) ", total);
  // Since first member may be skipped, do not print if empty. 
  if (ioarray_.empty())
    mprintf("\n");
  else
# else
  mprintf(" %zu members written) ", ioarray_.size());
# endif
    ioarray_.front()->Info();
  Traj().CommonInfo();
}
