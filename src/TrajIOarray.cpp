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
  replica_filenames_ = File::SearchForReplicas( tnameIn, debug_ );
  mprintf("\tFound %u replicas.\n", replica_filenames_.size());
  if (replica_filenames_.empty()) return 1;
  return 0;
}

/** Add lowest replica file name and names from comma-separated list. */
int TrajIOarray::AddReplicasFromArgs(FileName const& name0,
                                     std::string const& commaNames)
{
  if (name0.empty()) return 1;
  if (!File::Exists( name0 )) {
    File::ErrorMsg( name0.full() );
    return 1;
  }
  replica_filenames_.push_back( name0 );
  ArgList remdtraj_list( commaNames, "," );
  for (ArgList::const_iterator fname = remdtraj_list.begin();
                               fname != remdtraj_list.end(); ++fname)
  {
    FileName trajFilename( *fname );
    if (!File::Exists( trajFilename )) {
      File::ErrorMsg( trajFilename.full() );
      return 1;
    }
    replica_filenames_.push_back( trajFilename );
  }
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
    // Coordinates must be present.
    if (!replica0->CoordInfo().HasCrd()) {
      mprinterr("Error: No coordinates present in trajectory '%s'\n", repfile->full());
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
  replica_filenames_ = File::SearchForReplicas( tnameIn, trajComm.Master(), ensComm.Rank(),
                                                ensComm.Size(), debug_ );
  mprintf("\tFound %u replicas.\n", replica_filenames_.size());
  if (replica_filenames_.empty()) return 1;
  return 0;
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
      File::ErrorMsg( replica_filenames_[ensComm.Rank()].full() );
      rprinterr("Error: File '%s' not accessible.\n", replica_filenames_[ensComm.Rank()].full());
      return 1;
    }
  }
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
  mprintf("\tReading '%s' as %s\n", repFname.full(), TrajectoryFile::FormatString(repformat));
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
  // Coordinates must be present.
  if (!replica0->CoordInfo().HasCrd()) {
    mprinterr("Error: No coordinates present in trajectory '%s'\n", repFname.full());
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
  // Check # frames in all files, use lowest.
  Parallel::World().AllReduce( &totalFrames, &nframes, 1, MPI_INT, MPI_MIN );
  if (totalFrames != nframes) {
    rprintf("Warning: Replica '%s' frames (%i) is > # frames in shortest replica.\n",
            repFname.full(), nframes);
    mprintf("Warning: Setting total # of frames to read from replica ensemble to %i\n",
            totalFrames);
  }
  if (trajComm.Master()) {
    static const int iSize = 6;
    static const char* iTitle[iSize] = {"box", "velocity", "temperature", "time", "force",
                                        "replica dimensions"};
    // Check coordinate info of all files               0    1    2     3     4      5
    std::vector<int> Info( iSize * ensComm.Size() ); // box, vel, temp, time, force, nRepDims
    int rank_info[iSize];
    rank_info[0] = (int)cInfo.TrajBox().Type();
    rank_info[1] = (int)cInfo.HasVel();
    rank_info[2] = (int)cInfo.HasTemp();
    rank_info[3] = (int)cInfo.HasTime();
    rank_info[4] = (int)cInfo.HasForce();
    rank_info[5] = cInfo.ReplicaDimensions().Ndims();
    ensComm.AllGather( rank_info, iSize, MPI_INT, &Info[0] );
    // TODO Should mismatches be errors instead?
    for (int midx = 0; midx != iSize; midx++) {
      for (int ridx = midx + iSize; ridx < (int)Info.size(); ridx += iSize) {
        if (Info[midx] != Info[ridx]) {
          rprintf("Warning: Replica %i %s info does not match first replica.\n",
                  ridx/iSize, iTitle[midx]);
        }
      }
    }
  }
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
