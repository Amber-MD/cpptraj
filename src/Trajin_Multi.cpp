#include <locale> // isdigit
#include "Trajin_Multi.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, convertToInteger
#ifdef MPI
#  include "MpiRoutines.h"
#endif


// CONSTRUCTOR
Trajin_Multi::Trajin_Multi() :
  remdtrajtemp_(0.0),
  Ndimensions_(0),
  remd_indices_(0),
  lowestRepnum_(0),
  isSeekable_(true),
  hasVelocity_(false),
  isEnsemble_(false),
  replicasAreOpen_(false),
  targetType_(NONE),
  frameidx_(0)
# ifdef MPI
  ,ensembleFrameNum_(0)
# endif
{}

// DESTRUCTOR
Trajin_Multi::~Trajin_Multi() {
  if (replicasAreOpen_) EndTraj();
  for (IOarrayType::iterator replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    delete *replica;
  if (remd_indices_!=0) delete[] remd_indices_;
  if (frameidx_ != 0) delete[] frameidx_;
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
    mprinterr("Error: Traj %s has no numerical extension, required for automatic\n",
              TrajFilename().base());
    mprinterr("Error: detection of replica trajectories. Expected filename format is\n");
    mprinterr("Error: <Prefix>.<#> (with optional compression extension, examples:\n");
    mprinterr("Error: Rep.traj.nc.000,  remd.x.01.gz etc.\n");
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
  if (fileExists(replica_filename.c_str())) {
    mprintf("Warning: RemdTraj: Replica# found lower than file specified with trajin!\n");
    mprintf("Warning:           (Found %s)\n",replica_filename.c_str());
    mprintf("Warning:           trajin <file> remdtraj requires lowest # replica!\n");
  }

  // Add lowest filename, search for and add all replicas higher than it.
  ReplicaNames.push_back( TrajFilename().Full() );
  int current_repnum = lowestRepnum_;
  bool search_for_files = true;
  while (search_for_files) {
    ++current_repnum;
    replica_filename = Prefix + "." + integerToString(current_repnum, ExtWidth) + CompressExt;
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
  // Set base trajectory filename
  SetTrajFileName( tnameIn, true );
  // Process REMD-specific arguments
  if (argIn->Contains("remdtrajidx")) {
    // Looking for specific indices
    ArgList indicesArg(argIn->GetStringKey("remdtrajidx"), ",");
    if (indicesArg.empty()) {
      mprinterr("Error: remdtrajidx expects comma-separated list of target indices (e.g. 1,0,1)\n");
      return 1;
    }
    for (ArgList::const_iterator arg = indicesArg.begin(); 
                                 arg != indicesArg.end(); ++arg)
      remdtrajidx_.push_back( convertToInteger( *arg ) );
    targetType_ = INDICES;
  } else if (argIn->Contains("remdtrajtemp")) {
    // Looking for target temperature
    remdtrajtemp_ = argIn->getKeyDouble("remdtrajtemp",0.0);
    targetType_ = TEMP;
  }
  // If the command was ensemble, target args are not valid
  isEnsemble_ = false;
  bool no_sort = false;
  if ( argIn->CommandIs("ensemble") ){
    no_sort = argIn->hasKey("nosort");
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
  bool repBoxInfo = false;
  Ndimensions_ = -1;
  isSeekable_ = true;
  hasVelocity_ = false;
  for (NameListType::iterator repfile = replica_filenames_.begin();
                              repfile != replica_filenames_.end(); ++repfile)
  {
    // Detect format
    TrajectoryIO* replica0 = DetectFormat( *repfile );
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
      // Check how many frames will actually be read
      if (setupFrameInfo() == 0) return 1;
      // If lowest has box coords, check type against parm box info 
      Box parmBox = tparmIn->ParmBox();
      if (CheckBoxInfo(tparmIn->c_str(), parmBox, replica0->TrajBox())) return 1;
      tparmIn->SetBox( parmBox );
      // If lowest rep has box info, all others must have box info.
      repBoxInfo = replica0->HasBox();
      // If lowest rep has velocity, all others must have velocity.
      hasVelocity_ = replica0->HasV();
      // If lowest rep has dimensions, all others must have same dimensions
      Ndimensions_ = replica0->NreplicaDimensions();
      // Check that replica dimension valid for desired indices.
      if (targetType_ == INDICES && Ndimensions_ != (int)remdtrajidx_.size())
      {
        mprinterr("Error: RemdTraj: Replica # of dim (%i) not equal to target # dim (%zu)\n",
                  Ndimensions_, remdtrajidx_.size());
        return 1;
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
      if ( replica0->HasBox() != repBoxInfo ) {
        mprinterr("Error: RemdTraj: Replica %s box info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // Check velocity info against lowest rep.
      if ( replica0->HasV() != hasVelocity_ ) {
        mprinterr("Error: RemdTraj: Replica %s velocity info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
      // Check # dimensions against lowest rep
      if ( replica0->NreplicaDimensions() != Ndimensions_ ) {
        mprinterr("Error: RemdTraj: Replica %s dimension info does not match first replica.\n",
                  (*repfile).c_str());
        return 1;
      }
    }
    // All must be seekable or none will be
    if (isSeekable_ && !replica0->IsSeekable())
      isSeekable_ = false;
    // Check for temperature information. Not needed if not sorting.
    if ( !replica0->HasT() && !no_sort) {
      mprinterr("Error: RemdTraj: Replica %s does not have temperature info.\n",
                (*repfile).c_str());
      return 1;
    }
    lowestRep = false;
  }
  // Allocate space for number of dimensions
  if (Ndimensions_ > 0)
    remd_indices_ = new int[ Ndimensions_ ];
  // If targetType is currently NONE these will be processed as an ensemble. 
  // If dimensions are present index by replica indices, otherwise index
  // by temperature. If nosort was specified do not sort.
  if (isEnsemble_ && !no_sort) {
    if (Ndimensions_ > 0)
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
# ifdef MPI
  if (isEnsemble_) {
    // For ensemble, only open trajectory this thread will be dealing with
    //rprintf("Opening %s\n", replica_filenames_[worldrank].c_str()); // DEBUG
    if (REMDtraj_[worldrank]->openTrajin()) {
      rprinterr("Error: Trajin_Multi::BeginTraj: Could not open replica %s\n",
                replica_filenames_[worldrank].c_str());
      return 1;
    }
  } else {
# else
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
# endif
# ifdef MPI
  }
# endif
  // Set progress bar, start and offset.
  PrepareForRead( showProgress, isSeekable_ );
  replicasAreOpen_ = true;
  return 0;
}

// Trajin_Multi::EndTraj()
void Trajin_Multi::EndTraj() {
  if (replicasAreOpen_) {
#   ifdef MPI
    if (isEnsemble_)
      REMDtraj_[worldrank]->closeTraj();
    else
#   else
      for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
        (*replica)->closeTraj();
#   endif
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
    // Check if coords in frame are valid.
    if (frameIn.CheckCoordsInvalid())
      mprintf("Warning: Frame %i coords 1 & 2 overlap at origin; may be corrupt.\n",
              CurrentFrame()+1);
    tgtFrameFound = ProcessFrame();
  }

  return 1;
}

// Trajin_Multi::PrintInfo()
void Trajin_Multi::PrintInfo(int showExtended) {
  mprintf("REMD trajectories (%u total), lowest replica [%s]", REMDtraj_.size(),
          TrajFilename().base());
  if (showExtended == 1) PrintFrameInfo();
  mprintf("\n");
  if (debug_ > 0) {
    unsigned int repnum = 0;
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      mprintf("\t%u:[%s] ", repnum, replica_filenames_[repnum].c_str());
      ++repnum;
      (*replica)->Info();
      mprintf("\n");
    }
  }
  if (!isEnsemble_) {
    if (remdtrajidx_.empty())
      mprintf("\tLooking for frames at %.2lf K\n",remdtrajtemp_);
    else {
      mprintf("\tLooking for indices [");
      for (RemdIdxType::iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
        mprintf(" %i", *idx);
      mprintf(" ]\n");
    }
  } else {
    if ( targetType_ == INDICES )
      mprintf("\tProcessing ensemble using replica indices\n");
    else if ( targetType_ == TEMP )
      mprintf("\tProcessing ensemble using replica temperatures\n");
    else // NONE 
      mprintf("\tNot sorting ensemble.\n");
    if (debug_ > 0) EnsembleInfo();
  }
}

// -----------------------------------------------------------------------------
// Trajin_Multi::EnsembleInfo()
void Trajin_Multi::EnsembleInfo() const {
  if (targetType_ == TEMP) {
    mprintf("  Ensemble Temperature Map:\n");
    for (TmapType::const_iterator tmap = TemperatureMap_.begin();
                                  tmap != TemperatureMap_.end(); ++tmap)
      mprintf("\t%10.2f -> %i\n", (*tmap).first, (*tmap).second);
  } else if (targetType_ == INDICES) { // INDICES
    mprintf("  Ensemble Indices Map:\n");
    for (ImapType::const_iterator imap = IndicesMap_.begin();
                                  imap != IndicesMap_.end(); ++imap)
    {
      mprintf("\t{");
      for (std::vector<int>::const_iterator idx = (*imap).first.begin();
                                            idx != (*imap).first.end(); ++idx)
        mprintf(" %i", *idx);
      mprintf(" } -> %i\n", (*imap).second);
    }
  }
}

// Trajin_Multi::EnsembleSetup()
int Trajin_Multi::EnsembleSetup( FrameArray& f_ensemble ) {
  std::set<double> tList;
  std::set< std::vector<int> > iList;
  // Allocate space to hold position of each incoming frame in replica space.
  // TODO: When actually perfoming read in MPI will only need room for 1
  //frameidx_.resize( REMDtraj_.size() );
  if (frameidx_ != 0) delete[] frameidx_;
  frameidx_ = new int[ REMDtraj_.size() ];
  f_ensemble.resize( REMDtraj_.size() );
  f_ensemble.SetupFrames( TrajParm()->Atoms(), HasVelocity() );
  if (targetType_ == TEMP) {
    // Get a list of all temperature present in input REMD trajectories
    // by reading the first frame.
    // Assume that temperatures should be sorted lowest to highest.
    TemperatureMap_.clear();
    FrameArray::iterator frame = f_ensemble.begin();
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->openTrajin() ) return 1;
      if ( (*replica)->readFrame( CurrentFrame(), (*frame).xAddress(), (*frame).vAddress(),
                                  (*frame).bAddress(), (*frame).tAddress()) )
        return 1;
      (*replica)->closeTraj();
      std::pair<std::set<double>::iterator,bool> ret = tList.insert( (*frame).Temperature() );
      //std::pair<TmapType::iterator,bool> ret = 
      //  TemperatureMap_.insert(std::pair<double,int>((*frame).Temperature(),repnum++));
      if (!ret.second) {
        mprinterr("Error: Ensemble: Duplicate temperature detected (%.2f) in ensemble %s\n",
                  (*frame).Temperature(), TrajFilename().full());
        mprinterr("Info: If this is a H-REMD ensemble try the 'nosort' keyword.\n");
        return 1;
      }
    }
    // Temperatures are already sorted lowest to highest in set.
    int repnum = 0;
    for (std::set<double>::iterator temp0 = tList.begin(); temp0 != tList.end(); ++temp0)
      TemperatureMap_.insert(std::pair<double,int>(*temp0, repnum++)); 
  } else if (targetType_ == INDICES) {
    // Get a list of all indices present in input REMD trajectories
    // by reading in the first frame.
    IndicesMap_.clear();
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->openTrajin() ) return 1;
      if ( (*replica)->readIndices( CurrentFrame(), remd_indices_ ) ) return 1;
      (*replica)->closeTraj();
      std::pair<std::set< std::vector<int> >::iterator,bool> ret = 
        iList.insert( std::vector<int>( remd_indices_, remd_indices_ + Ndimensions_ ) );
      if ( !ret.second ) {
        mprinterr("Error: Ensemble: Duplicate indices detected in ensemble %s:",
                  TrajFilename().full());
        for (std::vector<int>::const_iterator idx = (*ret.first).begin(); 
                                              idx != (*ret.first).end(); ++idx)
          mprinterr(" %i", *idx);
        mprinterr("\n");
        return 1;
      }
    }
    int repnum = 0;
    for (std::set< std::vector<int> >::iterator idxs = iList.begin(); idxs != iList.end(); ++idxs)
      IndicesMap_.insert(std::pair< std::vector<int>, int >(*idxs, repnum++));
  } else { 
    // NONE, no sorting
#   ifdef MPI
    frameidx_[0] = worldrank;
#   else 
    for (int rnum = 0; rnum < (int) REMDtraj_.size(); ++rnum)
      frameidx_[rnum] = rnum;
#   endif
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
    //RemdIdxType::iterator fidx = frameidx_.begin();
    int* fidx = frameidx_;
    badEnsemble_ = false;
    // Read in all replicas
    //mprintf("DBG: Ensemble frame %i:",CurrentFrame()+1); // DEBUG
#   ifdef MPI
    // Read REMDtraj for this rank
    if ( REMDtraj_[worldrank]->readFrame( CurrentFrame(), (*frame).xAddress(), (*frame).vAddress(),
                                          (*frame).bAddress(), (*frame).tAddress()) )
      return 0;
#   else
    for (IOarrayType::iterator replica = REMDtraj_.begin(); replica!=REMDtraj_.end(); ++replica)
    {
      if ( (*replica)->readFrame( CurrentFrame(), (*frame).xAddress(), (*frame).vAddress(),
                                  (*frame).bAddress(), (*frame).tAddress()) )
        return 0;
#   endif
      if (targetType_ == TEMP) {
        TmapType::iterator tmap = TemperatureMap_.find( (*frame).Temperature() );
        if (tmap ==  TemperatureMap_.end())
          badEnsemble_ = true;
        else
          *fidx = (*tmap).second;
        //mprintf(" %.2f[%i]", (*frame).Temperature(), *fidx ); // DEBUG
      } else if (targetType_ == INDICES) {
#       ifdef MPI
        if ( REMDtraj_[worldrank]->readIndices( CurrentFrame(), remd_indices_ ) ) return 1;
#       else 
        if ( (*replica)->readIndices( CurrentFrame(), remd_indices_ ) ) return 1;
#       endif
        ImapType::iterator imap = IndicesMap_.find( 
          std::vector<int>( remd_indices_, remd_indices_ + Ndimensions_ ) );
        if (imap == IndicesMap_.end())
          badEnsemble_ = true;
        else
          *fidx = (*imap).second;
        // DEBUG
        //mprintf(" {");
        //for (int idx = 0; idx < Ndimensions_; ++idx)
        //  mprintf(" %i", remd_indices_[idx]);
        //mprintf(" }[%i]", *fidx);
      }
#     ifdef MPI
      // If calculated index is not worldrank, coords need to be sent to rank fidx.
      //rprintf("Index=%i\n", *fidx); // DEBUG
      ensembleFrameNum_ = 0;
      if (targetType_ != NONE) {
        // Each rank needs to know where to send its coords, and where to receive coords from.
        int my_idx = *fidx;
        if (parallel_allgather( &my_idx, 1, PARA_INT, frameidx_, 1, PARA_INT))
          rprinterr("Error: Gathering frame indices.\n");
        //mprintf("Frame %i Table:\n", CurrentFrame()); // DEBUG
        //for (unsigned int i = 0; i < REMDtraj_.size(); i++) // DEBUG
        //  mprintf("Rank %i has index %i\n", i, frameidx_[i]); // DEBUG
        // LOOP: one sendrecv at a time.
        for (int sendrank = 0; sendrank < (int)REMDtraj_.size(); sendrank++) {
          int recvrank = frameidx_[sendrank];
          if (sendrank != recvrank) {
            // TODO: Change Frame class so everything can be sent in one MPI call.
            if (sendrank == worldrank) {
              //rprintf("SENDING TO %i\n", recvrank); // DEBUG
              parallel_send( (*frame).xAddress(), (*frame).size(), PARA_DOUBLE, recvrank, 1212 );
              parallel_send( (*frame).bAddress(), 6, PARA_DOUBLE, recvrank, 1213 );
              parallel_send( (*frame).tAddress(), 1, PARA_DOUBLE, recvrank, 1214 );
              if (HasVelocity())
                parallel_send( (*frame).vAddress(), (*frame).size(), PARA_DOUBLE, recvrank, 1215 );
            } else if (recvrank == worldrank) {
              //rprintf("RECEIVING FROM %i\n", sendrank); // DEBUG
              parallel_recv( f_ensemble[1].xAddress(), (*frame).size(), PARA_DOUBLE, sendrank, 1212 );
              parallel_recv( f_ensemble[1].bAddress(), 6, PARA_DOUBLE, sendrank, 1213 );
              parallel_recv( f_ensemble[1].tAddress(), 1, PARA_DOUBLE, sendrank, 1214 );
              if (HasVelocity())
                parallel_recv( f_ensemble[1].vAddress(), (*frame).size(), PARA_DOUBLE, sendrank, 1215 );
              // Since a frame was received, indicate position 1 should be used
              ensembleFrameNum_ = 1; 
            }
          }
          //else rprintf("SEND RANK == RECV RANK, NO COMM\n"); // DEBUG
        } 
      }
      //rprintf("FRAME %i, FRAME RECEIVED= %i\n", CurrentFrame(), ensembleFrameNum_); // DEBUG 
#     else
      ++fidx;
      ++frame;
    }
#     endif
    //mprintf("\n"); // DEBUG
    tgtFrameFound = ProcessFrame();
  }
  return 1;
}
