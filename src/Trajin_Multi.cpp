#include <locale> // isdigit TODO: Use validInteger
#include "Trajin_Multi.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, convertToInteger
#include "DataFile.h" // For CRDIDX
#ifdef MPI
#  include "MpiRoutines.h"
#  ifdef TIMER
double Trajin_Multi::total_mpi_allgather_ = 0.0;
double Trajin_Multi::total_mpi_sendrecv_ = 0.0;
#  endif
#endif

// CONSTRUCTOR
Trajin_Multi::Trajin_Multi() :
  remdtrajtemp_(0.0),
  remdFrameFactor_(1.0),
  remdFrameOffset_(0),
  lowestRepnum_(0),
  replicasAreOpen_(false),
  targetType_(ReplicaInfo::NONE)
{}

// DESTRUCTOR
Trajin_Multi::~Trajin_Multi() {
  if (replicasAreOpen_) EndTraj();
  for (IOarrayType::iterator replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    delete *replica;
}

// Trajin_Multi::SearchForReplicas()
/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
Trajin_Multi::NameListType Trajin_Multi::SearchForReplicas() {
  NameListType ReplicaNames;
  std::string Prefix;
  std::string CompressExt;
  std::string ReplicaExt;
  std::locale loc;

  // STEP 1 - Get filename Prefix, Numerical extension, and optional
  //          compression extension.
  // Assume the extension of this trajectory is the number of the lowest 
  // replica, and that the other files are in sequence (e.g. rem.000, rem.001, 
  // rem.002 or rem.000.gz, rem.001.gz, rem.002.gz etc).
  if (debug_>1)
    mprintf("\tREMDTRAJ: FileName=[%s]\n",TrajFilename().full());
  if ( TrajFilename().Ext().empty() ) {
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n"
              "Error:   detection of replica trajectories. Expected filename format is\n"
              "Error:   <Prefix>.<#> (with optional compression extension), examples:\n"
              "Error:   Rep.traj.nc.000,  remd.x.01.gz etc.\n", TrajFilename().base());
    return ReplicaNames;
  }
  // Split off everything before replica extension
  size_t found = TrajFilename().Full().rfind( TrajFilename().Ext() );
  Prefix = TrajFilename().Full().substr(0, found); 
  ReplicaExt = TrajFilename().Ext(); // This should be the numeric extension
  // Remove leading '.'
  if (ReplicaExt[0] == '.') ReplicaExt.erase(0,1);
  CompressExt = TrajFilename().Compress();
  if (debug_>1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix.c_str(), ReplicaExt.c_str(), CompressExt.c_str());
  }

  // STEP 2 - Check that the numerical extension is valid.
  for (std::string::iterator schar = ReplicaExt.begin();
                             schar != ReplicaExt.end(); ++schar)
  {
    if (!isdigit(*schar,loc)) {
      mprinterr("Error: RemdTraj: Char [%c] in extension %s is not a digit!\n",
                *schar, ReplicaExt.c_str());
      return ReplicaNames;
    }
  }
  int ExtWidth = (int)ReplicaExt.size();
  if (debug_>1)
    mprintf("\tREMDTRAJ: Numerical Extension width=%i\n",ExtWidth);

  // STEP 3 - Store lowest replica number
  try { lowestRepnum_ = convertToInteger( ReplicaExt ); }
  catch ( ... ) {
    mprinterr("Error: RemdTraj: Could not convert lowest rep # %s to integer.\n",
              ReplicaExt.c_str());
    return ReplicaNames;
  }
  if (debug_>1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n",lowestRepnum_);

  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  std::string replica_filename = Prefix + "." + 
                                 integerToString(lowestRepnum_ - 1, ExtWidth) +
                                 CompressExt;
  if (fileExists(replica_filename)) {
    mprintf("Warning: Replica# found lower than file specified with trajin.\n"
            "Warning:   Found \"%s\"; 'trajin remdtraj' requires lowest # replica.\n",
            replica_filename.c_str());
  }

  // Add lowest filename, search for and add all replicas higher than it.
  ReplicaNames.push_back( TrajFilename().Full() );
  int current_repnum = lowestRepnum_;
  bool search_for_files = true;
  while (search_for_files) {
    ++current_repnum;
    replica_filename = Prefix + "." + integerToString(current_repnum, ExtWidth) + CompressExt;
    //mprintf("\t\tChecking for %s\n",replica_filename.c_str());
    if (fileExists(replica_filename))
      ReplicaNames.push_back( replica_filename );
    else
      search_for_files = false;
  }
  mprintf("\tFound %u replicas.\n",ReplicaNames.size());

  return ReplicaNames;
}

// Trajin_Multi::SetupTrajRead()
/** 'remdtraj' should have already been parsed out of the argIn list.
  */
int Trajin_Multi::SetupTrajRead(std::string const& tnameIn, ArgList& argIn, Topology *tparmIn)
{
  replica_filenames_.clear();
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Multi: No base filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Check for deprecated args
  if (argIn.hasKey("remdout")) {
    mprinterr("Error: 'remdout' is deprecated. To convert an entire replica ensemble the\n"
              "Error: correct usage is:\n"
              "Error:\t  ensemble <trajinfile> # (in place of 'trajin')\n"
              "Error:\t  trajout <trajoutfile> [<trajoutargs>]\n"
              "Error: Note that output trajectories will now have an integer appended to \n"
              "Error: them to indicate their place in the ensemble.\n");
    return 1;
  }
  // Check that base file exists
  if (!fileExists(tnameIn)) {
    mprinterr("Error: File %s does not exist.\n",tnameIn.c_str());
    return 1;
  }
  // Set base trajectory filename
  SetTrajFileName( tnameIn, true );
  // Process REMD-specific arguments
  if (argIn.Contains("remdtrajidx")) {
    // Looking for specific indices
    ArgList indicesArg(argIn.GetStringKey("remdtrajidx"), ",");
    if (indicesArg.empty()) {
      mprinterr("Error: remdtrajidx expects comma-separated list of target indices in each\n"
                "Error: dimension, '<dim1 idx>,<dim2 idx>,...,<dimN idx>'. Indices start\n"
                "Error: from 1.\n");
      return 1;
    }
    for (ArgList::const_iterator arg = indicesArg.begin(); 
                                 arg != indicesArg.end(); ++arg)
      remdtrajidx_.push_back( convertToInteger( *arg ) );
    targetType_ = ReplicaInfo::INDICES;
  } else if (argIn.Contains("remdtrajtemp")) {
    // Looking for target temperature
    remdtrajtemp_ = argIn.getKeyDouble("remdtrajtemp",0.0);
    targetType_ = ReplicaInfo::TEMP;
  }
  // Process any remlog keywords so they are not processed by SetupTrajIO
  std::string remlog_name = argIn.GetStringKey("remlog");
  double remlog_nstlim    = argIn.getKeyDouble("nstlim", 1.0);
  double remlog_ntwx      = argIn.getKeyDouble("ntwx",   1.0);
  // If the command was ensemble, target args are not valid
  bool no_sort = false;
  if ( IsEnsemble() ){
    no_sort = argIn.hasKey("nosort");
    if (targetType_ != ReplicaInfo::NONE || argIn.hasKey("remdtraj")) {
      mprintf("Warning: 'ensemble' does not use 'remdtraj', 'remdtrajidx' or 'remdtrajtemp'\n");
      targetType_ = ReplicaInfo::NONE;
    }
  } else {
    if (!remlog_name.empty()) {
      mprinterr("Error: 'remlog' is only for ensemble processing.\n");
      return 1;
    }
  }
  // CRDIDXARG: Parse out 'crdidx <indices list>' now so it is not processed
  //            by SetupTrajIO.
  ArgList crdidxarg;
  if (argIn.Contains("crdidx"))
    crdidxarg.SetList( "crdidx " + argIn.GetStringKey("crdidx"), " " );
  // Check if replica trajectories are explicitly listed
  ArgList remdtraj_list( argIn.GetStringKey("trajnames"), "," );
  if (remdtraj_list.Nargs()==0) {
    // Automatically scan for additional REMD traj files.
    replica_filenames_ = SearchForReplicas();
  } else {
    // Get filenames from args of remdtraj_list
    replica_filenames_.push_back( tnameIn );
    for (ArgList::const_iterator fname = remdtraj_list.begin();
                                 fname != remdtraj_list.end(); ++fname) 
       replica_filenames_.push_back( *fname );
  }
  
  // Loop over all filenames in replica_filenames 
  bool lowestRep = true;
  int Ndimensions = -1;
  cInfo_ = CoordinateInfo();
  TrajFormatType rep0format = TrajectoryFile::UNKNOWN_TRAJ;
  for (NameListType::iterator repfile = replica_filenames_.begin();
                              repfile != replica_filenames_.end(); ++repfile)
  {
    // Detect format
    TrajectoryIO* replica0 = DetectFormat( *repfile, rep0format );
    if ( replica0 == 0 ) {
      mprinterr("Error: RemdTraj: Could not set up replica file %s\n", (*repfile).c_str());
      return 1;
    }
    replica0->SetDebug( debug_ );
    // Pushing replica0 here allows the destructor to handle it on errors
    REMDtraj_.push_back( replica0 );
    if (lowestRep) {
      // Set up the lowest for reading and get the number of frames.
      if (SetupTrajIO( *repfile, *replica0, argIn )) return 1;
      cInfo_ = replica0->CoordInfo();
      Ndimensions = cInfo_.ReplicaDimensions().Ndims();
      // If lowest has box coords, check type against parm box info
      // If lowest rep has box info, all others must have box info.
      // If lowest rep has velocity, all others must have velocity.
      // If lowest rep has dimensions, all others must have same dimensions
      // Check that replica dimension valid for desired indices.
      if (targetType_ == ReplicaInfo::INDICES && Ndimensions != (int)remdtrajidx_.size())
      {
        mprinterr("Error: RemdTraj: Replica # of dim (%i) not equal to target # dim (%zu)\n",
                  Ndimensions, remdtrajidx_.size());
        return 1;
      }
      if (Ndimensions > 0) {
        mprintf("\tReplica dimensions:\n");
        for (int rd = 0; rd < Ndimensions; rd++)
          mprintf("\t\t%i: %s\n", rd+1, cInfo_.ReplicaDimensions().Description(rd)); 
      }
    } else {
      // Check total frames in this replica against lowest rep.
      int repframes = replica0->setupTrajin( *repfile, TrajParm() );
      if (repframes < 0 || repframes != TotalFrames()) {
        mprintf("Warning: RemdTraj: Replica %s frames (%i) does not match\n",
                 (*repfile).c_str(), repframes);
        mprintf("Warning:\t# frames in first replica (%i).\n",TotalFrames());
        if (repframes < 0) {
          mprinterr("Error: RemdTraj: Unknown # of frames in replica.\n");
          return 1;
        }
        if (repframes < TotalFrames()) {
          SetTotalFrames( repframes );
          mprintf("Warning: RemdTraj: Setting total # of frames to %i\n", TotalFrames());
        }
      }
      // Check box info against lowest rep.
      if ( replica0->CoordInfo().HasBox() != cInfo_.HasBox() ) {
        mprinterr("Error: RemdTraj: Replica %s box info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // TODO: Check specific box type
      // Check velocity info against lowest rep.
      if ( replica0->CoordInfo().HasVel() != cInfo_.HasVel() ) {
        mprinterr("Error: RemdTraj: Replica %s velocity info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // Check # dimensions and types against lowest rep
      if ( replica0->CoordInfo().ReplicaDimensions() != cInfo_.ReplicaDimensions() ) {
        mprinterr("Error: RemdTraj: Replica %s dimension info does not match first replica.\n",
                  (*repfile).c_str());
        ReplicaDimArray const& thisRepDims = replica0->CoordInfo().ReplicaDimensions();
        for (int rd = 0; rd < Ndimensions; rd++)
          mprinterr("\t\t%i: %s\n", rd+1, thisRepDims.Description(rd)); 
        return 1;
      }
      // If temperature/time info does not match set to false.
      if (cInfo_.HasTemp() != replica0->CoordInfo().HasTemp())
        cInfo_.SetTemperature( false );
      if (cInfo_.HasTime() != replica0->CoordInfo().HasTime())
        cInfo_.SetTime( false );
    }
    // Check for temperature information. Not needed if not sorting.
    if ( !replica0->CoordInfo().HasTemp() && !no_sort) {
      mprinterr("Error: RemdTraj: Replica %s does not have temperature info.\n",
                (*repfile).c_str());
      return 1;
    }
    lowestRep = false;
  }
  // Check how many frames will actually be read
  if (setupFrameInfo() == 0) return 1;
  // Unless nosort was specified, figure out target type if this will be 
  // processed as an ensemble.
  if (IsEnsemble() && !no_sort) {
    if ( !remlog_name.empty() ) {
      // Sort according to remlog data.
      DataFile remlogFile;
      DataSetList tempDSL;
      // CRDIDXARG: TODO: Come up with a way to do this that doesnt require ArgLists.
      if (remlogFile.ReadDataIn( remlog_name, crdidxarg, tempDSL ) ||
          tempDSL.empty())
      {
        mprinterr("Error: Could not read remlog data.\n");
        return 1;
      }
      if (remlogFile.Type() != DataFile::REMLOG) {
        mprinterr("Error: remlog: File was not of type remlog.\n");
        return 1;
      }
      if ( REMDtraj_.size() != tempDSL[0]->Size() ) {
        mprinterr("Error: ensemble size %zu does not match remlog ensemble size %zu\n",
                  REMDtraj_.size(), tempDSL[0]->Size());
        return 1;
      }
      remlogData_ = *((DataSet_RemLog*)tempDSL[0]); // FIXME: This feels clunky. Can we read direct?
      targetType_ = ReplicaInfo::CRDIDX;
      remdFrameFactor_ = remlog_ntwx / remlog_nstlim;
      mprintf("\t%g exchanges for every trajectory frame written.\n", remdFrameFactor_);
      if (remdFrameFactor_ > 1.0)
        remdFrameOffset_ = (int)remdFrameFactor_ - 1;
      else
        remdFrameOffset_ = 0;
      mprintf("\tTrajectory frame 1 corresponds to exchange %i\n", remdFrameOffset_ + 1);
      int expectedTrajFrames = (int)((double)TotalFrames() * remdFrameFactor_);
      if ( expectedTrajFrames != remlogData_.NumExchange() ) {
        mprinterr("Error: expected length of REMD ensemble %i does not match # exchanges in remlog %i.\n",
                  expectedTrajFrames, remlogData_.NumExchange());
        return 1;
      }
    } else {
      // If dimensions are present index by replica indices, otherwise index
      // by temperature. 
      if (Ndimensions > 0)
        targetType_ = ReplicaInfo::INDICES;
      else
        targetType_ = ReplicaInfo::TEMP;
    }
  }
  if (REMDtraj_.empty()) {
    mprinterr("Error: No replica trajectories set up.\n");
    return 1;
  }
  if (REMDtraj_.size() != replica_filenames_.size()) { // SANITY CHECK
    mprinterr("Error: Not all replica files were set up.\n");
    return 1;
  }
  // Update ensemble size
  cInfo_.SetEnsembleSize( (int)REMDtraj_.size() );
  if (debug_ > 0)
    Frame::PrintCoordInfo( TrajFilename().base(), TrajParm()->c_str(), cInfo_ );
  // FIXME: Should this ever be done here?
  TrajParm()->SetParmCoordInfo( cInfo_ ); 

  return 0;
}

// Trajin_Multi::BeginTraj()
int Trajin_Multi::BeginTraj(bool showProgress) {
# ifdef MPI
  if (IsEnsemble()) {
    // For ensemble, only open trajectory this thread will be dealing with
    //rprintf("Opening %s\n", replica_filenames_[worldrank].c_str()); // DEBUG
    if (REMDtraj_[worldrank]->openTrajin()) {
      rprinterr("Error: Trajin_Multi::BeginTraj: Could not open replica %s\n",
                replica_filenames_[worldrank].c_str());
      return 1;
    }
  } else {
# endif 
    // Open the trajectories
    mprintf("\tREMD: OPENING %zu REMD TRAJECTORIES\n", REMDtraj_.size());
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->openTrajin()) {
        mprinterr("Error: Trajin_Multi::BeginTraj: Could not open replica # %zu\n",
                  replica - REMDtraj_.begin() );
        return 1;
      }
    }
# ifdef MPI
  }
# endif
  // Set progress bar, start and offset.
  PrepareForRead( showProgress );
  replicasAreOpen_ = true;
  return 0;
}

// Trajin_Multi::EndTraj()
void Trajin_Multi::EndTraj() {
  if (replicasAreOpen_) {
#   ifdef MPI
    if (IsEnsemble()) {
      REMDtraj_[worldrank]->closeTraj();
#     ifdef TIMER
      total_mpi_allgather_ += mpi_allgather_timer_.Total();
      total_mpi_sendrecv_  += mpi_sendrecv_timer_.Total();
#     endif
    } else
#   endif 
      for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
        (*replica)->closeTraj();
    replicasAreOpen_ = false;
  }
}

#ifdef MPI
#ifdef TIMER
// Trajin_Multi::TimingData()
void Trajin_Multi::TimingData(double trajin_time) {
  if (total_mpi_allgather_ > 0.0 || total_mpi_sendrecv_ > 0.0) {
    double other_time = trajin_time - total_mpi_allgather_ - total_mpi_sendrecv_;
    rprintf("MPI_TIME:\tallgather: %.4f s (%.2f%%), sendrecv: %.4f s (%.2f%%), Other:  %.4f s (%.2f%%)\n",
            total_mpi_allgather_, (total_mpi_allgather_ / trajin_time)*100.0,
            total_mpi_sendrecv_,  (total_mpi_sendrecv_  / trajin_time)*100.0,
            other_time, (other_time / trajin_time)*100.0 );
  }
}
#endif
#endif

// Trajin_Multi::IsTarget()
/** Determine if given frame is target. ReadTrajFrame (i.e. non-ensemble) only.
   */
bool Trajin_Multi::IsTarget(Frame const& fIn) {
  if ( targetType_ == ReplicaInfo::TEMP ) {
    if ( fIn.Temperature() == remdtrajtemp_ ) return true;
  } else {
    Frame::RemdIdxType::const_iterator tgtIdx = fIn.RemdIndices().begin();
    for (RemdIdxType::const_iterator idx = remdtrajidx_.begin(); 
                                     idx != remdtrajidx_.end(); 
                                   ++idx, ++tgtIdx)
    {
      if ( *tgtIdx != *idx ) return false;
    }
    return true;
  }
  return false;
}

int Trajin_Multi::ReadTrajFrame( int currentFrame, Frame& frameIn ) {
  bool replicaFound = false;
  for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
  {
    // Locate the target temp/indices out of all the replicas
    if ( (*replica)->readFrame(currentFrame, frameIn))
      return 1;
    // Check if this is the target replica
    if ( IsTarget(frameIn) ) {
      replicaFound = true;
      break;
    }
  } // END loop over replicas
  if (!replicaFound) {
    mprinterr("Error: Target replica not found. Check that all replica trajectories\n"
              "Error:   were found and that the target temperature or indices are valid\n"
              "Error:   for this ensemble.\n");
    return 1; 
  }
  return 0;
}

// Trajin_Multi::PrintInfo()
void Trajin_Multi::PrintInfo(int showExtended) const {
  mprintf("REMD trajectories (%u total), lowest replica '%s'", REMDtraj_.size(),
          TrajFilename().base());
  if (showExtended == 1) PrintFrameInfo();
  mprintf("\n");
  if (debug_ > 0) {
    unsigned int repnum = 0;
    for (IOarrayType::const_iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      mprintf("\t%u:[%s] ", repnum, replica_filenames_[repnum].c_str());
      ++repnum;
      (*replica)->Info();
      mprintf("\n");
    }
  }
  if (!IsEnsemble()) {
    if (remdtrajidx_.empty())
      mprintf("\tLooking for frames at %.2lf K\n",remdtrajtemp_);
    else {
      mprintf("\tLooking for indices [");
      for (RemdIdxType::const_iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
        mprintf(" %i", *idx);
      mprintf(" ]\n");
    }
  } else {
    if ( targetType_ == ReplicaInfo::INDICES )
      mprintf("\tProcessing ensemble using replica indices\n");
    else if ( targetType_ == ReplicaInfo::TEMP )
      mprintf("\tProcessing ensemble using replica temperatures\n");
    else if ( targetType_ == ReplicaInfo::CRDIDX )
      mprintf("\tProcessing ensemble using remlog data, sorting by coordinate index.\n");
    else // NONE 
      mprintf("\tNot sorting ensemble.\n");
    if (debug_ > 0) EnsembleInfo();
  }
}

// -----------------------------------------------------------------------------
// Trajin_Multi::EnsembleInfo()
void Trajin_Multi::EnsembleInfo() const {
  if (targetType_ == ReplicaInfo::TEMP)
    PrintReplicaTmap( TemperatureMap_ );
  else if (targetType_ == ReplicaInfo::INDICES)
    PrintReplicaImap( IndicesMap_ );
  else if (targetType_ == ReplicaInfo::CRDIDX)
    mprintf("  Ensemble will be sorted by coordinate indices from remlog data.\n");
}

// Trajin_Multi::EnsembleSetup()
int Trajin_Multi::EnsembleSetup( FrameArray& f_ensemble, FramePtrArray& f_sorted ) {
  // Allocate space to hold position of each incoming frame in replica space.
# ifdef MPI
  // Only two frames needed; one for reading, one for receiving.
  f_sorted.resize( 2 );
  f_ensemble.resize( 2 );
  // This array will let each thread know who has what frame.
  frameidx_.resize( REMDtraj_.size() );
# else
  f_sorted.resize( REMDtraj_.size() );
  f_ensemble.resize( REMDtraj_.size() );
# endif
  f_ensemble.SetupFrames( TrajParm()->Atoms(), cInfo_ );
  // Get a list of all temperatures/indices.
  TemperatureMap_.ClearMap();
  IndicesMap_.ClearMap();
  if (targetType_ == ReplicaInfo::TEMP || targetType_ == ReplicaInfo::INDICES )
  {
#   ifdef MPI
    int err = 0;
    if (REMDtraj_[worldrank]->openTrajin()) {
      rprinterr("Error: Opening %s\n", TrajFilename().base());
      err = 1;
    }
    if (err == 0) {
      if (REMDtraj_[worldrank]->readFrame( Start(), f_ensemble[0] )) {
        rprinterr("Error: Reading %s\n", TrajFilename().base());
        err = 2;
      }
      REMDtraj_[worldrank]->closeTraj();
    }
    int total_error;
    parallel_allreduce( &total_error, &err, 1, PARA_INT, PARA_SUM );
    if (total_error != 0) {
      mprinterr("Error: Cannot setup ensemble trajectories.\n");
      return 1;
    }
#   else
    for (unsigned int member = 0; member != REMDtraj_.size(); member++) {
      if ( REMDtraj_[member]->openTrajin() ) return 1;
      if ( REMDtraj_[member]->readFrame( Start(), f_ensemble[member] ) ) return 1;
      REMDtraj_[member]->closeTraj();
    }
#   endif
    if (targetType_ == ReplicaInfo::TEMP) {
      TemperatureMap_ = SetReplicaTmap(REMDtraj_.size(), f_ensemble);
      if (TemperatureMap_.empty()) return 1;
    } else if (targetType_ == ReplicaInfo::INDICES) {
      IndicesMap_ = SetReplicaImap(REMDtraj_.size(),cInfo_.ReplicaDimensions().Ndims(),f_ensemble);
      if (IndicesMap_.empty()) return 1;
    }
  }  // Otherwise NONE, no sorting
  return 0;
}

// Trajin_Multi::ReadEnsemble()
int Trajin_Multi::ReadEnsemble( int currentFrame, FrameArray& f_ensemble, 
                                FramePtrArray& f_sorted )
{
  int fidx = 0;
  badEnsemble_ = false;
  // Read in all replicas
  //mprintf("DBG: Ensemble frame %i:",currentFrame+1); // DEBUG
# ifdef MPI
  int repIdx = worldrank; // for targetType==CRDIDX
  unsigned int member = 0;
  // Read REMDtraj for this rank
  if ( REMDtraj_[worldrank]->readFrame( currentFrame, f_ensemble[0]) )
    return 1;
# else
  int repIdx = 0; // for targetType==CRDIDX
  for (unsigned int member = 0; member != REMDtraj_.size(); ++member)
  {
    if ( REMDtraj_[member]->readFrame( currentFrame, f_ensemble[member]) )
      return 1;
# endif
    if (targetType_ == ReplicaInfo::TEMP)
      fidx = TemperatureMap_.FindIndex( f_ensemble[member].Temperature() );
    else if (targetType_ == ReplicaInfo::INDICES)
      fidx = IndicesMap_.FindIndex( f_ensemble[member].RemdIndices() );
    else if (targetType_ == ReplicaInfo::CRDIDX) {
      int currentRemExchange = (int)((double)currentFrame * remdFrameFactor_) + remdFrameOffset_;
      //mprintf("DEBUG:\tTrajFrame#=%i  RemdExch#=%i\n", currentFrame+1, currentRemExchange+1);
      fidx = remlogData_.RepFrame( currentRemExchange, repIdx++ ).CoordsIdx() - 1;
      //mprintf("DEBUG:\tFrame %i\tPosition %u is assigned index %i\n", currentFrame, member, fidx);
    }
#   ifndef MPI
    else // NONE
      fidx = member; // Not needed when MPI
    if (fidx == -1)
      badEnsemble_ = true;
    else
      f_sorted[fidx] = &f_ensemble[member];
  }
#   else 
  // If calculated index is not worldrank, coords need to be sent to rank fidx.
  //rprintf("Index=%i\n", fidx); // DEBUG
  int ensembleFrameNum = 0;
  if (targetType_ != ReplicaInfo::NONE) {
    // Each rank needs to know where to send its coords, and where to receive coords from.
#   ifdef TIMER
    mpi_allgather_timer_.Start();
#   endif
    if (parallel_allgather( &fidx, 1, PARA_INT, &frameidx_[0], 1, PARA_INT)) {
      rprinterr("Error: Gathering frame indices.\n");
      return 0; // TODO: Better parallel error check
    }
    for (unsigned int i = 0; i != REMDtraj_.size(); i++)
      if (frameidx_[i] == -1) {
        badEnsemble_ = true;
        break;
      }
#   ifdef TIMER
    mpi_allgather_timer_.Stop();
    mpi_sendrecv_timer_.Start();
#   endif
    // LOOP: one sendrecv at a time.
    if (!badEnsemble_) {
      for (int sendrank = 0; sendrank != (int)REMDtraj_.size(); sendrank++) {
        int recvrank = frameidx_[sendrank];
        if (sendrank != recvrank) {
          if (sendrank == worldrank)
            f_ensemble[0].SendFrame( recvrank ); 
          else if (recvrank == worldrank) {
            f_ensemble[1].RecvFrame( sendrank );
            // Since a frame was received, indicate position 1 should be used
            ensembleFrameNum = 1; 
          }
        }
        //else rprintf("SEND RANK == RECV RANK, NO COMM\n"); // DEBUG
      }
    }
#   ifdef TIMER
    mpi_sendrecv_timer_.Stop();
#   endif
  }
  f_sorted[0] = &f_ensemble[ensembleFrameNum];
  //rprintf("FRAME %i, FRAME RECEIVED= %i\n", currentFrame, ensembleFrameNum); // DEBUG 
#   endif
  return 0;
}

/** CRDIDXARG:
  * \return A string containing the coordinate indices (comma separated) of the
  *         final exchange in remlogData_.
  */
std::string Trajin_Multi::FinalCrdIndices() const {
  if (remlogData_.Empty()) return std::string();
  std::string arg("crdidx ");
  int finalExchg = remlogData_.NumExchange() - 1;
  for (unsigned int rep = 0; rep < remlogData_.Size(); rep++) {
    if (rep > 0) arg += ",";
    arg += ( integerToString( remlogData_.RepFrame(finalExchg, rep).CoordsIdx() ) );
  }
  return arg;
}
