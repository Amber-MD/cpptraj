#include "Trajout_Multi.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename
#ifdef MPI
# include "MpiRoutines.h"
#endif

// CONSTRUCTOR
Trajout_Multi::Trajout_Multi() : ensembleSize_(0) {} 

// DESTRUCTOR
Trajout_Multi::~Trajout_Multi() {
  EndTraj();
  Clear();
}

void Trajout_Multi::Clear() {
  for (IOarrayType::const_iterator tio = ioarray_.begin(); tio != ioarray_.end(); ++tio)
    delete *tio;
  ioarray_.clear();
}

// Trajout_Multi::InitTrajWrite()
/** Initialize each output trajectory with appropriate TrajectoryIO class
  * and process arguments.
  */
int Trajout_Multi::InitTrajWrite(std::string const& tnameIn, ArgList const& argIn, 
                                 Topology *tparmIn, TrajFormatType writeFormatIn)
{
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: InitTrajWrite: No filename given.\n");
    return 1;
  }
  ensembleSize_ = tparmIn->ParmCoordInfo().EnsembleSize();
  // Require that ensemble size is set.
  if (ensembleSize_ < 1) {
    mprinterr("Internal Error: Ensemble size too small for ensemble output.\n");
    return 1;
  }
  ArgList trajout_args = argIn;
  TrajectoryFile::TrajFormatType writeFormat = writeFormatIn;
  // Get onlymembers range
  Range members_to_write = MembersToWrite(trajout_args.GetStringKey("onlymembers"), ensembleSize_);
  if (members_to_write.Empty()) return 1;
  // DEBUG
  //std::string dbg_mtw = "MembersToWrite:";
  //for (Range::const_iterator r = members_to_write.begin(); r != members_to_write.end(); r++)
  //  dbg_mtw += (" " + integerToString(*r));
  //rprintf("DEBUG: %s\n", dbg_mtw.c_str());
  // Process common args
  if (CommonTrajoutSetup(tnameIn, trajout_args, tparmIn, writeFormat))
    return 1;
  Clear();
  // Set up ensemble file names.
  fileNames_.clear();
# ifdef MPI
  // In MPI each thread writes a single member.
  if (members_to_write.InRange( worldrank ))
    fileNames_.push_back( NumberFilename(TrajFilename().Full(), worldrank) );
  else
    rprintf("Warning: Skipping member '%s'\n", 
            NumberFilename(TrajFilename().Full(), worldrank).c_str());
# else
  // In serial single process writes each member.
  // Create a map: tIndex[ pos ] = <ioarray_index>
  tIndex_.clear();
  tIndex_.reserve( ensembleSize_ );
  int ioidx = 0;
  for (int num = 0; num < ensembleSize_; num++) {
    if (members_to_write.InRange( num )) {
      fileNames_.push_back( NumberFilename(TrajFilename().Full(), num) );
      tIndex_.push_back( ioidx++ );
    } else {
      mprintf("Warning: Skipping member '%s'\n",
              NumberFilename(TrajFilename().Full(), num).c_str());
      tIndex_.push_back( -1 );
    }
  }
# endif
  // Set up write format for each file. 
  typedef std::vector<TrajectoryFile::TrajFormatType> FmtArray;
  FmtArray fileFormats(fileNames_.size(), writeFormat);
  // If appending, all files must exist and must have same format.
  if (TrajoutAppend()) {
    for (unsigned int m = 0; m != fileNames_.size(); ++m)
      CheckAppendFormat( fileNames_[m], fileFormats[m] );
  }
  // Set up TrajectoryIO for each member.
  for (unsigned int m = 0; m != fileNames_.size(); ++m) {
    rprintf("\tWriting ensemble '%s' as %s\n", fileNames_[m].c_str(),
            TrajectoryFile::FormatString(fileFormats[m]));
    TrajectoryIO* tio = AllocTrajIO( fileFormats[m] );
    if (tio == 0) return 1;
    ioarray_.push_back( tio );
    ioarray_.back()->SetDebug( debug_ );
    // Set specified title - will not set if empty 
    //if (!TrajoutTitle().empty())
    //  trajio_->SetTitle( TrajoutTitle() + "." + integerToString( num ) );
    ioarray_.back()->SetTitle( TrajoutTitle() );
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

// Trajout_Multi::EndTraj()
void Trajout_Multi::EndTraj() {
  if (TrajIsOpen()) {
    for (IOarrayType::const_iterator tio = ioarray_.begin(); tio != ioarray_.end(); ++tio)
      (*tio)->closeTraj();
    SetTrajIsOpen( false );
  }
}

/** Perform any topology-related setup for this trajectory if given Topology
  * matches what trajectory was initialized with; the topology may have
  * been modified (e.g. by a 'strip' command) since the output trajectory was
  * initialized.
  */
int Trajout_Multi::SetupTrajWrite(Topology* tparmIn) {
  // First frame setup
  if (!TrajIsOpen()) {
    for (unsigned int m = 0; m != ioarray_.size(); ++m) {
      if (FirstFrameSetup(fileNames_[m], ioarray_[m], tparmIn)) return 1;
    }
  }
  return 0;
}

// Trajout_Multi::WriteEnsemble()
/** Write given array of frames if trajectory is open (initialzed and set-up).
  */ 
int Trajout_Multi::WriteEnsemble(int set, FramePtrArray const& Farray)
{
  // Check that set should be written
  if (CheckFrameRange(set)) return 0;
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

// Trajout_Multi::PrintInfo()
void Trajout_Multi::PrintInfo(int showExtended) const {
  mprintf("  '%s.X' (Ensemble,", TrajFilename().base());
# ifdef MPI
  // Since not every thread may be writing if 'onlymembers' specified,
  // determine total number being written.
  int mysize = (int)ioarray_.size();
  int total;
  parallel_reduce(&total, &mysize, 1, PARA_INT, PARA_SUM);
  mprintf(" %i members written) ", total);
  // Since first member may be skipped, do not print if empty. 
  if (ioarray_.empty())
    mprintf("\n");
  else
# else
  mprintf(" %zu members written) ", ioarray_.size());
# endif
  CommonInfo( ioarray_.front() );
}
