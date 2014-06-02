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
  // Require that ensemble size is set.
  if (ensembleSize_ < 1) {
    mprinterr("Internal Error: Ensemble size has not been set.\n");
    return 1;
  }
  ArgList trajout_args = argIn;
  TrajectoryFile::TrajFormatType writeFormat = writeFormatIn;
  // Process common args
  if (CommonTrajoutSetup(tnameIn, trajout_args, tparmIn, writeFormat))
    return 1;
  Clear();
  // Set up ensemble file names.
  fileNames_.clear();
# ifdef MPI
  // In MPI each thread writes a single member.
  fileNames_.push_back( NumberFilename(TrajFilename().Full(), worldrank) );
# else
  // In serial single process writes each member.
  for (int num = 0; num < total; num++)
    fileNames_.push_back( NumberFilename(TrajFilename().Full(), num) );
# endif
  // Set up write format for each file. 
  typedef std::vector<TrajectoryFile::TrajFormatType> FmtArray;
  FmtArray fileFormats(fileNames_.size(), writeFormat);
  // If appending, all files must exist and must have same format.
  if (TrajoutAppend()) {
    FmtArray::iterator wfmt = fileFormats.begin(); 
    for (Sarray::const_iterator fname = fileNames_.begin();
                                fname != fileNames_.end(); 
                              ++fname, ++wfmt)
    {
      CheckAppendFormat( *fname, *wfmt );
    }
  }
  // FIXME: Should report heterogeneous ensembles
  mprintf("\tWriting ensemble '%s.X' as %s\n", TrajFilename().full(),
          TrajectoryFile::FormatString(writeFormat));
  // Set up TrajectoryIO for each member.
  FmtArray::const_iterator wfmt;
  for (Sarray::const_iterator fname = fileNames_.begin();
                              fname != fileNames_.end(); 
                            ++fname, ++wfmt)
  {
    TrajectoryIO* tio = AllocTrajIO( *wfmt );
    if (tio == 0) return 1;
    ioarray_.push_back( tio );
    ioarray_.back()->SetDebug( debug_ );
    // Set specified title - will not set if empty 
    //if (!TrajoutTitle().empty())
    //  trajio_->SetTitle( TrajoutTitle() + "." + integerToString( num ) );
    ioarray_.back()->SetTitle( TrajoutTitle() );
    // Process any write arguments specific to certain formats not related
    // to parm file. Options related to parm file are handled on the first
    // write in WriteFrame.
    ArgList rep_args = trajout_args;
    if (ioarray_.back()->processWriteArgs( rep_args )) {
      mprinterr("Error: trajout %s: Could not process arguments.\n",fname->c_str());
      return 1;
    }
  }
  // Write is set up for topology only when first frame written.
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

// Trajout_Multi::WriteEnsemble()
/** Write given array of frames to their correct ensemble position. If this is
  * the first frame being written, this routine is where the output trajectories
  * will actually be set up for the associated topology file since the topology 
  * may have been modified (e.g. by a 'strip' command) since the output 
  * trajectory was initialized (modified topologies will still have the 
  * same Pindex).
  */ 
int Trajout_Multi::WriteEnsemble(int set, Topology *tparmIn, FrameArray const& Farray,
                                 Frame::RemdIdxType const& FrameIdx)
{
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->Pindex() != TrajParm()->Pindex()) return 0;

  // First frame setup
  if (!TrajIsOpen()) {
    Sarray::const_iterator fname = fileNames_.begin();
    for (IOarrayType::const_iterator tio = ioarray_.begin(); 
                                     tio != ioarray_.end(); 
                                   ++tio, ++fname)
      if (FirstFrameSetup(*fname, *tio, tparmIn)) return 1;
  }

  // Check that set should be written
  if (CheckFrameRange(set)) return 0;

  // Write
  for (int member = 0; member != ensembleSize_; member++) {
    mprintf("DEBUG: member %i to position %i\n", member, FrameIdx[member]);
    if (ioarray_[ FrameIdx[member] ]->writeFrame(set, Farray[member])) return 1;
  }
  
  return 0;
}

// Trajout_Multi::PrintInfo()
void Trajout_Multi::PrintInfo(int showExtended) const {
  mprintf("  (Ensemble)");
  CommonInfo( ioarray_.front() );
}
