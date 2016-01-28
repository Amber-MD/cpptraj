#include "TrajIOarray.h"
#include "StringRoutines.h" // integerToString, validInteger
#include "CpptrajStdio.h"
#include "TrajectoryFile.h"

TrajIOarray::~TrajIOarray() { ClearIOarray(); }

const char* TrajIOarray::DEPRECATED_remdout = 
  "Error: 'remdout' is deprecated. To convert an entire replica ensemble the\n"
  "Error: correct usage is:\n"
  "Error:\t  ensemble <trajinfile> # (in place of 'trajin')\n"
  "Error:\t  trajout <trajoutfile> [<trajoutargs>]\n"
  "Error: Note that output trajectories will now have an integer appended to \n"
  "Error: them to indicate their place in the ensemble.\n";

void TrajIOarray::ClearIOarray() {
  for (IOarrayType::const_iterator rep = IOarray_.begin(); rep != IOarray_.end(); ++rep)
    if (*rep != 0) delete *rep; // TODO: Close?
  IOarray_.clear();
  replica_filenames_.clear();
}

int TrajIOarray::SetupReplicaFilenames(FileName const& tnameIn, ArgList& argIn) {
  std::string trajnames = argIn.GetStringKey("trajnames");
  if (!trajnames.empty())
    return AddReplicasFromArgs( tnameIn, trajnames );
  return SearchForReplicas( tnameIn );
}

/** Add lowest replica file name and names from comma-separated list. */
int TrajIOarray::AddReplicasFromArgs(FileName const& name0,
                                     std::string const& commaNames)
{
  if (name0.empty()) return 1;
  if (!File::Exists( name0 )) {
    mprinterr("Error: File '%s' does not exist.\n", name0.full());
    return 1;
  }
  replica_filenames_.push_back( name0 );
  ArgList remdtraj_list( commaNames, "," );
  for (ArgList::const_iterator fname = remdtraj_list.begin();
                               fname != remdtraj_list.end(); ++fname)
  {
    FileName trajFilename( *fname );
    if (!File::Exists( trajFilename )) {
      mprinterr("Error: File '%s' does not exist.\n", trajFilename.full());
      return 1;
    }
    replica_filenames_.push_back( trajFilename );
  }
  return 0;
}

/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
int TrajIOarray::SearchForReplicas(FileName const& fname) {
  RepName repName(fname, debug_);
  if (repName.Error()) return 1;
  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  if (File::Exists( repName.RepFilename( -1 ) )) {
    mprintf("Warning: Replica# found lower than file specified with trajin.\n"
            "Warning:   Found \"%s\"; 'trajin remdtraj' requires lowest # replica.\n",
            repName.RepFilename( -1 ).full());
  }
  // Add lowest replica filename, search for and add all replicas higher than it.
  replica_filenames_.push_back( fname );
  int rep_offset = 0;
  bool search_for_files = true;
  FileName trajFilename;
  while (search_for_files) {
    ++rep_offset;
    trajFilename = repName.RepFilename( rep_offset );
    //mprintf("\t\tChecking for %s\n", trajFilename.full());
    if (File::Exists( trajFilename ))
      replica_filenames_.push_back( trajFilename );
    else
      search_for_files = false;
  }
  mprintf("\tFound %u replicas.\n", replica_filenames_.size());

  return 0;
}

/** Loop over all filenames in replica_filenames, set up TrajectoryIO. */
int TrajIOarray::SetupIOarray(ArgList& argIn, TrajFrameCounter& counter,
                              CoordinateInfo& cInfo, Topology* trajParm)
{
  // Sanity check
  if (!IOarray_.empty()) {
    mprinterr("Internal Error: SetupIOarray() has been called twice.\n");
    return 1;
  }
  // Save arguments that have not been processed so they can be passed
  // to each replica in turn. Only the lowest replica will use argIn.
  ArgList argCopy( argIn );
  bool lowestRep = true;
  int rep0Frames = TrajectoryIO::TRAJIN_UNK;  // Total frames in replica 0
  int totalFrames = TrajectoryIO::TRAJIN_UNK; // Total # frames to be read from ensemble
  TrajectoryFile::TrajFormatType lastRepFmt = TrajectoryFile::UNKNOWN_TRAJ;
  // Right now enforce that all replicas have the same metadata as lowest
  // replica, e.g. if replica 0 has temperature, replica 1 does too etc.
  for (File::NameArray::const_iterator repfile = replica_filenames_.begin();
                                       repfile != replica_filenames_.end(); ++repfile)
  {
    // Detect format
    TrajectoryFile::TrajFormatType repformat = TrajectoryFile::UNKNOWN_TRAJ;
    TrajectoryIO* replica0 = TrajectoryFile::DetectFormat( *repfile, repformat );
    if ( replica0 == 0 ) {
      mprinterr("Error: Could not set up replica file %s\n", repfile->full());
      return 1;
    }
    if (repformat != lastRepFmt)
      mprintf("\tReading '%s' as %s\n", repfile->full(), TrajectoryFile::FormatString(repformat));
    lastRepFmt = repformat;
    replica0->SetDebug( debug_ );
    // Pushing replica0 here allows the destructor to handle it on errors
    IOarray_.push_back( replica0 );
    // Process format-specific read args. Do not exit on error in case
    // replicas have different formats supporting different args.
    if (lowestRep) {
      replica0->processReadArgs( argIn );
    } else {
      ArgList argtmp( argCopy );
      replica0->processReadArgs( argtmp );
    }
    // Set up replica for reading and get the number of frames.
    int nframes = replica0->setupTrajin( *repfile, trajParm );
    if (nframes == TrajectoryIO::TRAJIN_ERR) {
      mprinterr("Error: Could not set up %s for reading.\n", repfile->full());
      return 1;
    }
    // TODO: Do not allow unknown number of frames?
    if (lowestRep) {
      cInfo = replica0->CoordInfo();
      rep0Frames = nframes;
      totalFrames = nframes;
      if (cInfo.ReplicaDimensions().Ndims() > 0) {
        mprintf("\tReplica dimensions:\n");
        for (int rd = 0; rd < cInfo.ReplicaDimensions().Ndims(); rd++)
          mprintf("\t\t%i: %s\n", rd+1, cInfo.ReplicaDimensions().Description(rd));
      }
    } else {
      // Check total frames in this replica against lowest rep.
      if (nframes != rep0Frames)
        mprintf("Warning: Replica %s frames (%i) does not match # frames in first replica (%i).\n",
                repfile->base(), nframes, rep0Frames);
      //if (repframes < 0) {
      //  mprinterr("Error: RemdTraj: Unknown # of frames in replica.\n");
      //  return 1;
      //}
      if (nframes < totalFrames) {
        totalFrames = nframes; 
        mprintf("Warning: Setting total # of frames to read from replica ensemble to %i\n",
                totalFrames);
      }
      // Check box info against lowest rep.
      if ( replica0->CoordInfo().HasBox() != cInfo.HasBox() ) {
        mprinterr("Error: Replica %s box info does not match first replica.\n",
                  repfile->full());
        return 1;
      }
      // TODO: Check specific box type
      // Check velocity info against lowest rep.
      if ( replica0->CoordInfo().HasVel() != cInfo.HasVel() ) {
        mprinterr("Error: Replica %s velocity info does not match first replica.\n",
                  repfile->full());
        return 1;
      }
      // Check # dimensions and types against lowest rep
      if ( replica0->CoordInfo().ReplicaDimensions() != cInfo.ReplicaDimensions() ) {
        mprinterr("Error: Replica %s dimension info does not match first replica.\n",
                  repfile->full());
        ReplicaDimArray const& thisRepDims = replica0->CoordInfo().ReplicaDimensions();
        for (int rd = 0; rd < thisRepDims.Ndims(); rd++)
          mprinterr("\t\t%i: %s\n", rd+1, thisRepDims.Description(rd)); 
        return 1;
      }
      // If temperature/time info does not match set to false.
      if (cInfo.HasTemp() != replica0->CoordInfo().HasTemp())
        cInfo.SetTemperature( false );
      if (cInfo.HasTime() != replica0->CoordInfo().HasTime())
        cInfo.SetTime( false );
    }
    lowestRep = false;
  }
  // Check how many frames will actually be read
  if (counter.CheckFrameArgs( totalFrames, argIn )) return 1;
  // Check for errors.
  if (IOarray_.empty()) {
    mprinterr("Error: No replica trajectories set up.\n");
    return 1;
  }
  if (IOarray_.size() != replica_filenames_.size()) { // SANITY CHECK
    mprinterr("Error: Not all replica files were set up.\n");
    return 1;
  }
  // Update ensemble size
  cInfo.SetEnsembleSize( (int)IOarray_.size() );
  if (debug_ > 0)
    cInfo.PrintCoordInfo( replica_filenames_[0].full(), trajParm->c_str() );

  return 0;
}

void TrajIOarray::PrintIOinfo() const {
  for (unsigned int rn = 0; rn != IOarray_.size(); ++rn) {
    mprintf("\t%u:[%s] ", rn, replica_filenames_[rn].base());
    if (IOarray_[rn] != 0) IOarray_[rn]->Info();
    mprintf("\n");
  }
}
#ifdef MPI
// ----- PARALLEL ROUTINES -----------------------------------------------------
/** Setup replica filenames in parallel. */
int TrajIOarray::SetupReplicaFilenames(FileName const& tnameIn, ArgList& argIn,
                                       Parallel::Comm const& ensComm,
                                       Parallel::Comm const& trajComm)
{
  std::string trajnames = argIn.GetStringKey("trajnames");
  if (!trajnames.empty())
    return AddReplicasFromArgs( tnameIn, trajnames, ensComm, trajComm );
  return SearchForReplicas( tnameIn, ensComm, trajComm );
}

/** Each rank checks that specified file is present. */
int TrajIOarray::AddReplicasFromArgs(FileName const& name0,
                                     std::string const& commaNames,
                                     Parallel::Comm const& ensComm,
                                     Parallel::Comm const& trajComm)
{
  // First set up filename array on all ranks.
  if (name0.empty()) return 1;
  replica_filenames_.push_back( name0 );
  ArgList remdtraj_list( commaNames, "," );
  for (ArgList::const_iterator fname = remdtraj_list.begin();
                               fname != remdtraj_list.end(); ++fname)
    replica_filenames_.push_back( FileName( *fname ) );
  if (ensComm.Size() != (int)replica_filenames_.size())
    return 1;
  else if (trajComm.Master()) { // Only traj comm master checks file
    if (!File::Exists( replica_filenames_[ ensComm.Rank() ])) {
      rprinterr("Error: File '%s' does not exist.\n", replica_filenames_[ensComm.Rank()].full());
      return 1;
    }
  }
  return 0;
}

/** Each rank searches for replica based on lowest replica number. */
int TrajIOarray::SearchForReplicas(FileName const& fname, Parallel::Comm const& ensComm,
                                   Parallel::Comm const& trajComm)
{
  RepName repName(fname, debug_);
  if (repName.Error()) return 1;
  // TODO check for lower replica number?
  FileName replicaFilename = repName.RepFilename( ensComm.Rank() );
  // Only traj comm masters actually check for files.
  if (trajComm.Master()) {
    if (!File::Exists( replicaFilename )) return 1;
  }
  // At this point each rank has found its replica. Populate filename array.
  for (int offset = 0; offset < ensComm.Size(); ++offset)
    replica_filenames_.push_back( repName.RepFilename( offset ) );
  mprintf("\tFound %zu replicas.\n", replica_filenames_.size());
  return 0;
}

/** Each rank only sets up file that it will process. */
int TrajIOarray::SetupIOarray(ArgList& argIn, TrajFrameCounter& counter,
                              CoordinateInfo& cInfo, Topology* trajParm,
                              Parallel::Comm const& ensComm, Parallel::Comm const& trajComm)
{
  // Sanity check
  if (!IOarray_.empty()) {
    mprinterr("Internal Error: SetupIOarray() has been called twice.\n");
    return 1;
  }
  // Detect format
  FileName const& repFname = replica_filenames_[ensComm.Rank()];
  TrajectoryFile::TrajFormatType repformat = TrajectoryFile::UNKNOWN_TRAJ;
  TrajectoryIO* replica0 = TrajectoryFile::DetectFormat( repFname, repformat );
  if ( replica0 == 0 ) {
    mprinterr("Error: Could not set up replica file %s\n", repFname.full());
    return 1;
  }
  rprintf("\tReading '%s' as %s\n", repFname.full(), TrajectoryFile::FormatString(repformat));
  replica0->SetDebug( debug_ );
  // Construct the IOarray_ with blanks for all except this rank.
  for (int member = 0; member != ensComm.Size(); member++)
    if (member == ensComm.Rank())
      IOarray_.push_back( replica0 );
    else
      IOarray_.push_back( 0 );
  // Process format-specific read args.
  replica0->processReadArgs( argIn );
  // Set up replica for reading and get # frames
  int nframes = replica0->setupTrajin( repFname, trajParm );
  if (nframes == TrajectoryIO::TRAJIN_ERR) {
    mprinterr("Error: Could not set up %s for reading.\n", repFname.full());
    return 1;
  }
  // Set coordinate info
  cInfo = replica0->CoordInfo();
  int totalFrames = nframes;
  if (cInfo.ReplicaDimensions().Ndims() > 0) { // TODO put in common routine
    mprintf("\tReplica dimensions:\n");
    for (int rd = 0; rd < cInfo.ReplicaDimensions().Ndims(); rd++)
      mprintf("\t\t%i: %s\n", rd+1, cInfo.ReplicaDimensions().Description(rd));
  }
  // TODO: Check coordinate info of all files
  // TODO: Put code below into a common routine with serial version
  // Check how many frames will actually be read
  if (counter.CheckFrameArgs( totalFrames, argIn )) return 1;
  // SANITY CHECK
  if (IOarray_.size() != replica_filenames_.size()) {
    mprinterr("Error: Not all replica files were set up.\n");
    return 1;
  }
  // Update ensemble size
  cInfo.SetEnsembleSize( (int)IOarray_.size() );
  if (debug_ > 0)
    cInfo.PrintCoordInfo( repFname.full(), trajParm->c_str() );

  return 0;
}
#endif

// ----- RepName Class ---------------------------------------------------------
// RepName CONSTRUCTOR
TrajIOarray::RepName::RepName(FileName const& fname, int debugIn) {
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: FileName=[%s]\n", fname.full());
  if ( fname.Ext().empty() ) {
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n"
              "Error:   detection of replica trajectories. Expected filename format is\n"
              "Error:   <Prefix>.<#> (with optional compression extension), examples:\n"
              "Error:   Rep.traj.nc.000,  remd.x.01.gz etc.\n", fname.base());
    return;
  }
  // Split off everything before replica extension
  size_t found = fname.Full().rfind( fname.Ext() );
  Prefix_.assign( fname.Full().substr(0, found) );
  ReplicaExt_.assign( fname.Ext() ); // This should be the numeric extension
  // Remove leading '.'
  if (ReplicaExt_[0] == '.') ReplicaExt_.erase(0,1);
  CompressExt_.assign( fname.Compress() );
  if (debugIn > 1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix_.c_str(), ReplicaExt_.c_str(), CompressExt_.c_str());
  }
  // Check that the numerical extension is valid.
  if ( !validInteger(ReplicaExt_) ) {
    mprinterr("Error: Replica extension [%s] is not an integer.\n", ReplicaExt_.c_str());
    Prefix_.clear(); // Empty Prefix_ indicates error.
    return;
  }
  ExtWidth_ = (int)ReplicaExt_.size();
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: Numerical Extension width=%i\n", ExtWidth_);
  // Store lowest replica number
  lowestRepnum_ = convertToInteger( ReplicaExt_ );
  // TODO: Do not allow negative replica numbers?
  if (debugIn > 1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n", lowestRepnum_);
}

/** \return Replica file name for given offset from lowest replica number. */
FileName TrajIOarray::RepName::RepFilename(int offset) const {
  FileName trajFilename;
  trajFilename.SetFileName_NoExpansion( Prefix_ + "." +
                                        integerToString(lowestRepnum_ + offset, ExtWidth_) +
                                        CompressExt_ );
  return trajFilename;
}
