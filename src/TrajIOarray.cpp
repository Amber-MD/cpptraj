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
    delete *rep; // TODO: Close?
  IOarray_.clear();
  replica_filenames_.clear();
}

int TrajIOarray::SetupReplicaFilenames(FileName const& tnameIn, ArgList& argIn) {
  std::string trajnames = argIn.GetStringKey("trajnames");
  if (!trajnames.empty()) {
    if (AddReplicasFromArgs( tnameIn, trajnames )) return 1;
  } else {
    if (SearchForReplicas( tnameIn )) return 1;
  }
  return 0;
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

#ifdef MPI
/** Setup replica filenames in parallel. */
int TrajIOarray::SetupReplicaFilenames(FileName const& tnameIn, ArgList& argIn,
                                       Parallel::Comm const& ensComm,
                                       Parallel::Comm const& trajComm)
{
  std::string trajnames = argIn.GetStringKey("trajnames");
  if (trajComm.Master()) {
    if (!trajnames.empty()) {
      if (AddReplicasFromArgs( tnameIn, trajnames, ensComm )) return 1;
    } else {
      if (SearchForReplicas( tnameIn, ensComm )) return 1;
    }
  }
  // Make sure all trajectory ranks have file names.
  if (trajComm.Size() > 1) {
    int eSize;              // # ensemble members
    size_t bufSize;         // Total size of the character buffer.
    std::vector<char> cbuf; // Filename character buffer
    if (trajComm.Master()) {
      // Broadcast total number of replicas
      eSize = (int)replica_filenames_.size();
      trajComm.MasterBcast( &eSize, 1, MPI_INT );
      // Determine/broadcast the buffer size needed for all replica file names + null chars
      bufSize = 0;
      for (File::NameArray::const_iterator fn = replica_filenames_.begin();
                                           fn != replica_filenames_.end(); ++fn)
      {
        bufSize += fn->Full().size();
        bufSize++; // null
      }
      cbuf.resize( bufSize );
      trajComm.MasterBcast( &bufSize, 1, MPI_UNSIGNED_LONG_LONG );
      // Put all file names into buffer separated by nulls and broadcast.
      std::vector<char>::iterator ptr = cbuf.begin();
      for (int idx = 0; idx != eSize; idx++) {
        std::copy( replica_filenames_[idx].Full().begin(),
                   replica_filenames_[idx].Full().end(), ptr );
        ptr += replica_filenames_[idx].Full().size();
        *(ptr++) = '\0';
      }
      trajComm.MasterBcast( &cbuf[0], bufSize, MPI_CHAR );
    } else {
      // Receive total number of replicas
      trajComm.MasterBcast( &eSize, 1, MPI_INT );
      replica_filenames_.reserve( eSize );
      // Receive buffer size needed for all replica file names + null chars
      trajComm.MasterBcast( &bufSize, 1, MPI_UNSIGNED_LONG_LONG );
      cbuf.resize( bufSize );
      // Receive all file names in a buffer and put into replica_filenames_
      trajComm.MasterBcast( &cbuf[0], bufSize, MPI_CHAR );
      std::vector<char>::const_iterator beg = cbuf.begin();
      for (int idx = 0; idx != eSize; idx++) {
        // Seek to null char
        std::vector<char>::const_iterator end = beg;
        while ( *end != '\0' ) ++end;
        replica_filenames_.push_back( std::string(beg, end) );
        beg = end + 1;
      }
    }
  }
  return 0;
}

/** Each rank checks that specified file is present. */
int TrajIOarray::AddReplicasFromArgs(FileName const& name0,
                                     std::string const& commaNames,
                                     Parallel::Comm const& commIn)
{
  // First set up filename array on all ranks.
  if (name0.empty()) return 1;
  replica_filenames_.push_back( name0 );
  ArgList remdtraj_list( commaNames, "," );
  for (ArgList::const_iterator fname = remdtraj_list.begin();
                               fname != remdtraj_list.end(); ++fname)
    replica_filenames_.push_back( FileName( *fname ) );
  int err = 0;
  if (commIn.Size() != (int)replica_filenames_.size())
    err = 1;
  else if (!File::Exists( replica_filenames_[ commIn.Rank() ])) {
    rprinterr("Error: File '%s' does not exist.\n", replica_filenames_[commIn.Rank()].full());
    err = 1;
  }
  if (commIn.CheckError( err )) return 1;

  return 0;
}

/** Each rank searches for replica based on lowest replica number. */
int TrajIOarray::SearchForReplicas(FileName const& fname, Parallel::Comm const& commIn) {
  RepName lowestRepName(fname, debug_);
  if (lowestRepName.Error()) return 1;
  // TODO check for higher replica number?
  FileName replicaFilename = lowestRepName.RepFilename( commIn.Rank() +
                                                        lowestRepName.LowestRepNum() );
  int err = 0;
  if (!File::Exists( replicaFilename )) err = 1;
  if (commIn.CheckError( err )) return 1;
  rprintf("DEBUG: Found '%s'\n", replicaFilename.full());
  // At this point each rank has found its replica. Populate filename array.
  for (int repnum = lowestRepName.LowestRepNum();
           repnum < lowestRepName.LowestRepNum() + commIn.Size(); ++repnum)
    replica_filenames_.push_back( lowestRepName.RepFilename( repnum ) );
  mprintf("\tFound %zu replicas.\n", replica_filenames_.size());
  return 0;
}
#endif

/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
int TrajIOarray::SearchForReplicas(FileName const& fname) {
  // STEP 1 - Get filename Prefix, Numerical extension, and optional
  //          compression extension.
  // Assume the extension of this trajectory is the number of the lowest 
  // replica, and that the other files are in sequence (e.g. rem.000, rem.001, 
  // rem.002 or rem.000.gz, rem.001.gz, rem.002.gz etc).
  if (debug_>1)
    mprintf("\tREMDTRAJ: FileName=[%s]\n",fname.full());
  if ( fname.Ext().empty() ) {
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n"
              "Error:   detection of replica trajectories. Expected filename format is\n"
              "Error:   <Prefix>.<#> (with optional compression extension), examples:\n"
              "Error:   Rep.traj.nc.000,  remd.x.01.gz etc.\n", fname.base());
    return 1;
  }
  // Split off everything before replica extension
  size_t found = fname.Full().rfind( fname.Ext() );
  std::string Prefix = fname.Full().substr(0, found); 
  std::string ReplicaExt = fname.Ext(); // This should be the numeric extension
  // Remove leading '.'
  if (ReplicaExt[0] == '.') ReplicaExt.erase(0,1);
  std::string CompressExt = fname.Compress();
  if (debug_>1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix.c_str(), ReplicaExt.c_str(), CompressExt.c_str());
  }

  // STEP 2 - Check that the numerical extension is valid.
  if ( !validInteger(ReplicaExt) ) {
    mprinterr("Error: Replica extension [%s] is not an integer.\n", ReplicaExt.c_str());
    return 1;
  }
  int ExtWidth = (int)ReplicaExt.size();
  if (debug_>1)
    mprintf("\tREMDTRAJ: Numerical Extension width=%i\n",ExtWidth);

  // STEP 3 - Store lowest replica number
  int lowestRepnum = convertToInteger( ReplicaExt );
  // TODO: Do not allow negative replica numbers?
  if (debug_>1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n",lowestRepnum);
  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  std::string replica_filename = Prefix + "." + 
                                 integerToString(lowestRepnum - 1, ExtWidth) +
                                 CompressExt;
  if (File::Exists(replica_filename)) {
    mprintf("Warning: Replica# found lower than file specified with trajin.\n"
            "Warning:   Found \"%s\"; 'trajin remdtraj' requires lowest # replica.\n",
            replica_filename.c_str());
  }

  // SETP 4 - Add lowest filename, search for and add all replicas higher than it.
  replica_filenames_.push_back( fname );
  int current_repnum = lowestRepnum;
  bool search_for_files = true;
  FileName trajFilename;
  while (search_for_files) {
    ++current_repnum;
    trajFilename.SetFileName_NoExpansion( Prefix + "." +
                                          integerToString(current_repnum, ExtWidth) +
                                          CompressExt );
    //mprintf("\t\tChecking for %s\n",replica_filename.c_str());
    if (File::Exists(trajFilename))
      replica_filenames_.push_back( trajFilename );
    else
      search_for_files = false;
  }
  mprintf("\tFound %u replicas.\n",replica_filenames_.size());

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
    IOarray_[rn]->Info();
    mprintf("\n");
  }
}
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

/** \return Replica file name for given replica number. */
FileName TrajIOarray::RepName::RepFilename(int repnum) const {
  FileName trajFilename;
  trajFilename.SetFileName_NoExpansion( Prefix_ + "." + integerToString(repnum, ExtWidth_) +
                                        CompressExt_ );
  return trajFilename;
}
