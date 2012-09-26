#include <locale> // isdigit
#include <sstream>
#include "Trajin_Multi.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Trajin_Multi::Trajin_Multi() :
  remdtrajtemp_(0.0),
  remd_indices_(0),
  lowestRepnum_(0),
  isSeekable_(true),
  hasVelocity_(false),
  isEnsemble_(false),
  replicasAreOpen_(false),
  targetType_(NONE)
{}

// DESTRUCTOR
Trajin_Multi::~Trajin_Multi() {
  if (replicasAreOpen_) EndTraj();
  for (IOarrayType::iterator replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    delete *replica;
}

// PrintNoExtError()
static void PrintNoExtError(const char* nameIn) {
  mprinterr("Error: Traj %s has no numerical extension, required for automatic\n", nameIn);
  mprinterr("       detection of replica trajectories. Expected filename format is\n");
  mprinterr("       <Prefix>.<#> (with optional compression extension, examples:\n");
  mprinterr("       Rep.traj.nc.000,  remd.x.01.gz etc.\n");
}

// Trajin_Multi::SearchForReplicas()
/** Assuming lowest replica filename has been set, search for all other 
  * replica names assuming a naming scheme of '<PREFIX>.<EXT>[.<CEXT>]', 
  * where <EXT> is a numerical extension and <CEXT> is an optional 
  * compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
Trajin_Multi::NameListType Trajin_Multi::SearchForReplicas(bool isCompressed) {
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
    mprintf("\tREMDTRAJ: FileName=[%s]\n",FullTrajStr());
  size_t found = FullTrajName().find_last_of(".");
  // NOTE: Eventually allow just numbers, no prefix? 000 002 etc
  if (found == std::string::npos) {
    PrintNoExtError(BaseTrajName().c_str());
    return ReplicaNames;
  }
  // If file is not compressed, this should be the numeric extension
  if (!isCompressed) {
    Prefix = FullTrajName().substr(0, found);
    ReplicaExt = FullTrajName().substr(found+1);
  } else {
  // If file is compressed, this should be the compression extension
  // include the dot.
    CompressExt = FullTrajName().substr(found);
    // Get prefix and numerical extension, without the Compression extension
    std::string compressPrefix = FullTrajName().substr(0,found);
    found = compressPrefix.find_last_of(".");
    if (found == std::string::npos) {
      PrintNoExtError(BaseTrajName().c_str());
      return ReplicaNames;
    }
    Prefix = compressPrefix.substr(0,found);
    ReplicaExt = compressPrefix.substr(found+1);
  }
  if (debug_>1) {
    mprintf("\tREMDTRAJ: Prefix=[%s], #Ext=[%s], CompressExt=[%s]\n",
            Prefix.c_str(), ReplicaExt.c_str(), CompressExt.c_str());
  }

  // STEP 2 - Check that the numerical extension is valid.
  for (std::string::iterator schar = ReplicaExt.begin();
                             schar != ReplicaExt.end();
                             schar++)
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
  std::istringstream iss( ReplicaExt );
  if ( !(iss >> lowestRepnum_) ) {
    mprinterr("Error: RemdTraj: Could not convert lowest rep # %s to integer.\n",
              ReplicaExt.c_str());
    return ReplicaNames;
  }
  if (debug_>1)
    mprintf("\tREMDTRAJ: index of first replica = %i\n",lowestRepnum_);

  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  std::ostringstream replica_num;
  replica_num.fill('0');
  replica_num.width( ExtWidth );
  replica_num << std::right << (lowestRepnum_ - 1);
  std::string replica_filename = Prefix + "." + replica_num.str() + CompressExt;
  if (fileExists(replica_filename.c_str())) {
    mprintf("    Warning: RemdTraj: Replica# found lower than file specified with trajin!\n");
    mprintf("             (Found %s)\n",replica_filename.c_str());
    mprintf("             trajin <file> remdtraj requires lowest # replica!\n");
  }

  // Add lowest filename, search for and add all replicas higher than it.
  ReplicaNames.push_back( FullTrajName() );
  int current_repnum = lowestRepnum_;
  bool search_for_files = true;
  while (search_for_files) {
    ++current_repnum;
    replica_num.str("");
    replica_num.fill('0');
    replica_num.width( ExtWidth );
    replica_num << std::right << current_repnum;
    replica_filename = Prefix + "." + replica_num.str() + CompressExt;
    //mprintf("\t\tChecking for %s\n",replica_filename.c_str());
    if (fileExists(replica_filename.c_str()))
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
int Trajin_Multi::SetupTrajRead(std::string const& tnameIn, ArgList *argIn, Topology *tparmIn)
{
  replica_filenames_.clear();
  // Require a base filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Multi: No base filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Parallel Trajectory processing currenly requires some args to set up
  if (argIn == 0) {
    mprinterr("Internal Error: Trajin_Multi: No arguments given.\n");
    return 1;
  }
  // Check for deprecated args
  if (argIn->hasKey("remdout")) {
    mprinterr("Error: 'remdout' is deprecated. To convert an entire replica ensemble\n");
    mprinterr("Error: the correct usage is:\n");
    mprinterr("Error:\t  ensemble <trajinfile> # (in place of 'trajin')\n");
    mprinterr("Error:\t  trajout <trajoutfile> [<trajoutargs>]\n");
    mprinterr("Error: Note that output trajectories will now have an integer\n");
    mprinterr("Error: appended to them to indicate their place in the ensemble.\n");
    return 1;
  }
  // Check that base file exists
  if (!fileExists(tnameIn.c_str())) {
    mprinterr("Error: File %s does not exist.\n",tnameIn.c_str());
    return 1;
  }
  // Set up CpptrajFile
  CpptrajFile baseFile;
  if (baseFile.SetupRead( tnameIn, debug_ )) return 1;
  // Set file name and base file name
  SetFileNames( tnameIn, baseFile.BaseFileName() );
  // Process REMD-specific arguments
  if (argIn->Contains("remdtrajidx")) {
    // Looking for specific indices
    ArgList indicesArg(argIn->GetStringKey("remdtrajidx"), ",");
    if (indicesArg.empty()) {
      mprinterr("Error: remdtrajidx expects comma-separated list of target indices (e.g. 1,0,1)\n");
      return 1;
    }
    for (int argidx = 0; argidx < indicesArg.Nargs(); ++argidx)
       remdtrajidx_.push_back( indicesArg.ArgToInteger(argidx) );
    targetType_ = INDICES;
    remd_indices_ = new int[ remdtrajidx_.size() ];
  } else if (argIn->Contains("remdtrajtemp")) {
    // Looking for target temperature
    remdtrajtemp_ = argIn->getKeyDouble("remdtrajtemp",0.0);
    targetType_ = TEMP;
  }
  // If the command was ensemble, target args are not valid
  isEnsemble_ = false;
  if ( argIn->CommandIs("ensemble") ){
    isEnsemble_ = true;
    if (targetType_ != NONE || argIn->hasKey("remdtraj")) {
      mprintf("Warning: 'ensemble' does not use 'remdtraj', 'remdtrajidx' or 'remdtrajtemp'\n");
      targetType_ = NONE;
    }
  }

  // Check if replica trajectories are explicitly listed
  ArgList remdtraj_list( argIn->GetStringKey("trajnames"), "," );
  if (remdtraj_list.Nargs()==0) {
    // Automatically scan for additional REMD traj files.
    replica_filenames_ = SearchForReplicas(baseFile.IsCompressed());
  } else {
    // Get filenames from args of remdtraj_list
    replica_filenames_.push_back( tnameIn );
    for (ArgList::const_iterator fname = remdtraj_list.begin();
                                 fname != remdtraj_list.end(); ++fname) 
       replica_filenames_.push_back( *fname );
  }

  // Loop over all filenames in replica_filenames 
  bool lowestRep = true;
  bool repBoxInfo = false;
  int Ndimensions = -1;
  isSeekable_ = true;
  hasVelocity_ = false;
  for (NameListType::iterator repfile = replica_filenames_.begin();
                              repfile != replica_filenames_.end(); ++repfile)
  {
    // Set file info
    if (!lowestRep)
      baseFile.SetupRead( *repfile, debug_ );
    // Detect format
    TrajectoryIO* replica0 = DetectFormat( baseFile );
    if ( replica0 == 0 ) {
      mprinterr("Error: RemdTraj: Could not set up replica file %s\n", (*repfile).c_str());
      return 1;
    }
    // Pushing replica0 here allows the destructor to handle it on errors
    REMDtraj_.push_back( replica0 );
    if (lowestRep) {
      // Set up the lowest for reading and get the number of frames.
      if (SetupTrajIO( replica0, argIn )) return 1;
      // Check how many frames will actually be read
      if (setupFrameInfo() == 0) return 1;
      // If lowest has box coords, set box type from box Angles
      if (replica0->CheckBoxInfo( TrajParm() ) != 0) {
        mprinterr("Error in replica %s box info.\n", (*repfile).c_str());
        return 1;
      }
      // If lowest rep has box info, all others must have box info.
      repBoxInfo = replica0->HasBox();
      // If lowest rep has velocity, all others must have velocity.
      hasVelocity_ = replica0->HasVelocity();
      // If lowest rep has dimensions, all others must have same dimensions
      Ndimensions = replica0->NreplicaDimensions();
      // Check that replica dimension valid for desired indices.
      if (targetType_ == INDICES && Ndimensions != (int)remdtrajidx_.size())
      {
        mprinterr("Error: RemdTraj: Replica # of dim (%i) not equal to target # dim (%zu)\n",
                  Ndimensions, remdtrajidx_.size());
        return 1;
      }
    } else {
      // Check total frames in this replica against lowest rep.
      int repframes = replica0->setupTrajin( TrajParm() );
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
      if ( replica0->HasBox() != repBoxInfo ) {
        mprinterr("Error: RemdTraj: Replica %s box info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // Check velocity info against lowest rep.
      if ( replica0->HasVelocity() != hasVelocity_ ) {
        mprinterr("Error: RemdTraj: Replica %s velocity info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // Check # dimensions against lowest rep
      if ( replica0->NreplicaDimensions() != Ndimensions ) {
        mprinterr("Error: RemdTraj: Replica %s dimension info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
    }
    // All must be seekable or none will be
    if (isSeekable_ && !replica0->Seekable())
      isSeekable_ = false;
    // Check for temperature information
    if ( !replica0->HasTemperature()) {
      mprinterr("Error: RemdTraj: Replica %s does not have temperature info.\n",
                (*repfile).c_str());
      return 1;
    }
    lowestRep = false;
  }
  // If targetType is currently NONE these will be processed as an ensemble. 
  // If dimensions are present index by replica indices, otherwise index
  // by temperature.
  if (isEnsemble_) {
    if (Ndimensions > 0)
      targetType_ = INDICES;
    else
      targetType_ = TEMP;
  }
  if (REMDtraj_.empty()) {
    mprinterr("Error: No replica trajectories set up.\n");
    return 1;
  }
  if (REMDtraj_.size() != replica_filenames_.size()) { // SANITY CHECK
    mprinterr("Error: Not all replica files were set up.\n");
    return 1;
  }
  if (!isSeekable_) {
    mprinterr("Error: Currently all replica trajectories must be seekable.\n");
    return 1;
  }
  
  return 0;
}

// Trajin_Multi::BeginTraj()
int Trajin_Multi::BeginTraj(bool showProgress) {
  // Open the trajectories
  mprintf("\tREMD: OPENING REMD TRAJECTORIES\n");
  for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
  {
    if ( (*replica)->openTraj()) {
      mprinterr("Error: Trajin_Multi::BeginTraj: Could not open replica # %zu\n",
                replica - REMDtraj_.begin() );
      return 1;
    }
  }
  // Set progress bar, start and offset.
  PrepareForRead( showProgress, isSeekable_ );
  replicasAreOpen_ = true;
  return 0;
}

// Trajin_Multi::EndTraj()
void Trajin_Multi::EndTraj() {
  if (replicasAreOpen_) {
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
      (*replica)->closeTraj();
    replicasAreOpen_ = false;
  }
}

// Trajin_Multi::IsTarget()
bool Trajin_Multi::IsTarget(double tempIn) {
  if ( targetType_ == TEMP ) {
    if ( tempIn == remdtrajtemp_ ) return true;
  } else {
    const int* tgtIdx = remd_indices_;
    for (RemdIdxType::iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
    {
      if ( *tgtIdx != *idx ) return false;
      ++tgtIdx;
    }
    return true;
  }
  return false;
}

// Trajin_Multi::GetNextFrame()
int Trajin_Multi::GetNextFrame( Frame& frameIn ) {
  // If the current frame is out of range, exit
  if ( CheckFinished() ) return 0;

  bool tgtFrameFound = false;
  while ( !tgtFrameFound ) {
    bool replicaFound = false;
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      // Locate the target temp/indices out of all the replicas
      if ( (*replica)->readFrame(CurrentFrame(), frameIn.xAddress(), frameIn.vAddress(),
                                 frameIn.bAddress(), frameIn.tAddress()))
        return 0;
      (*replica)->readIndices(CurrentFrame(),remd_indices_);
      // Check if this is the target replica
      if ( IsTarget(frameIn.Temperature()) ) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        //mprintf("REMDTRAJ: Replica %zu matches!\n", reptrajin - REMDtraj_.begin());
        // TODO: I think if !isSeekable this will break the read since some
        //       trajectories will not have readFrame called and so will be
        //       behind those that did.
        replicaFound = true;
        break;
      }
    } // END loop over replicas
    if (!replicaFound) {
      mprinterr("Error: Target replica not found. Check that all replica trajectories\n");
      mprinterr("Error: were found and that the target temperature or indices are valid\n");
      mprinterr("Error: for this ensemble.\n");
      return 0; 
    }
    tgtFrameFound = ProcessFrame();
  }

  return 1;
}

// Trajin_Multi::PrintInfo()
void Trajin_Multi::PrintInfo(int showExtended) {
  mprintf("  REMD trajectories (%u total), lowest replica [%s]", REMDtraj_.size(),
          BaseTrajStr());
  if (showExtended == 1) PrintFrameInfo();
  mprintf("\n");
  if (showExtended == 1) {
    unsigned int repnum = 0;
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      mprintf("\t%u:[%s] ", repnum, replica_filenames_[repnum].c_str());
      ++repnum;
      (*replica)->info();
      mprintf("\n");
    }
  }
  if (!isEnsemble_) {
    if (remdtrajidx_.empty())
      mprintf("\tLooking for frames at %.2lf K",remdtrajtemp_);
    else {
      mprintf("\tLooking for indices [");
      for (RemdIdxType::iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
        mprintf(" %i", *idx);
      mprintf(" ]");
    }
  } else {
    mprintf("\tProcessing ensemble using");
    if ( targetType_ == INDICES )
      mprintf(" replica indices\n");
    else {
      mprintf(" replica temperatures\n");
      mprintf("\tEnsemble Temperature Map:\n");
      for (TmapType::iterator tmap = TemperatureMap_.begin();
                                          tmap != TemperatureMap_.end(); ++tmap)
        mprintf("\t\t%10.2f -> %i\n", (*tmap).first, (*tmap).second);
    }
  }
}

// -----------------------------------------------------------------------------
// Trajin_Multi::EnsembleSetup()
int Trajin_Multi::EnsembleSetup( FrameArray& f_ensemble ) {
  // Allocate space to hold position of each incoming frame in replica space.
  frameidx_.resize( REMDtraj_.size() );
  f_ensemble.resize( REMDtraj_.size() );
  f_ensemble.SetupFrames( TrajParm()->Atoms(), HasVelocity() );
  if (targetType_ == TEMP) {
    // Get a list of all temperature present in input REMD trajectories
    // by reading the first frame.
    TemperatureMap_.clear();
    FrameArray::iterator frame = f_ensemble.begin();
    int repnum = 0;
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->openTraj() ) return 1;
      if ( (*replica)->readFrame( CurrentFrame(), (*frame).xAddress(), (*frame).vAddress(),
                                  (*frame).bAddress(), (*frame).tAddress()) )
        return 1;
      (*replica)->closeTraj();
      std::pair<TmapType::iterator,bool> ret = 
        TemperatureMap_.insert(std::pair<double,int>((*frame).Temperature(),repnum++));
      if (!ret.second) {
        mprinterr("Error: Ensemble: Duplicate temperature detected (%.2f)\n",
                  (*frame).Temperature());
        return 1;
      }
    } 
  } else if (targetType_ == INDICES) {
    return 1;
  } 
  return 0;
}

// Trajin_Multi::GetNextEnsemble()
int Trajin_Multi::GetNextEnsemble( FrameArray& f_ensemble ) {
  // If the current frame is out of range, exit
  if ( CheckFinished() ) return 0;
  bool tgtFrameFound = false;
  while ( !tgtFrameFound ) {
    FrameArray::iterator frame = f_ensemble.begin();
    RemdIdxType::iterator fidx = frameidx_.begin();
    // Read in all replicas
    //mprintf("DBG: Ensemble frame %i:",CurrentFrame()+1); // DEBUG
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->readFrame( CurrentFrame(), (*frame).xAddress(), (*frame).vAddress(),
                                  (*frame).bAddress(), (*frame).tAddress()) )
        return 0;
      // TODO: Indices read
      TmapType::iterator tmap = TemperatureMap_.find( (*frame).Temperature() );
      *fidx = (*tmap).second;
      //mprintf(" %.2f[%i]", (*frame).Temperature(), *fidx ); // DEBUG
      ++fidx;
      ++frame;
    }
    //mprintf("\n"); // DEBUG
    tgtFrameFound = ProcessFrame();
  }
  return 1;
}


