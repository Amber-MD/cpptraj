#include "Trajin_Multi.h"
#include "StringRoutines.h" // convertToInteger
#include "CpptrajStdio.h"

Trajin_Multi::Trajin_Multi() :
  targetType_(ReplicaInfo::NONE),
  remdtrajtemp_(0.0)
{}

/** Expect lowest replica file name to be passed in. 'remdtraj' should have
  * already been parsed out of input arguments.
  */
int Trajin_Multi::SetupTrajRead(FileName const& tnameIn, ArgList& argIn, Topology *tparmIn)
{
  // Set file name and topology pointer.
  if (SetTraj().SetNameAndParm(tnameIn, tparmIn)) return 1;
  REMDtraj_.ClearIOarray();
  // Check for deprecated args
  if (argIn.hasKey("remdout")) {
    mprinterr("%s", TrajIOarray::DEPRECATED_remdout);
    return 1;
  }
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
    for (ArgList::const_iterator arg = indicesArg.begin(); // TODO: validInteger? 
                                 arg != indicesArg.end(); ++arg)
      remdtrajidx_.push_back( convertToInteger( *arg ) );
    targetType_ = ReplicaInfo::INDICES;
  } else if (argIn.Contains("remdtrajtemp")) {
    // Looking for target temperature
    remdtrajtemp_ = argIn.getKeyDouble("remdtrajtemp",0.0);
    targetType_ = ReplicaInfo::TEMP;
  }
  // Set up replica file names.
  if (REMDtraj_.SetupReplicaFilenames( tnameIn, argIn )) return 1;

  // Set up TrajectoryIO classes for all file names.
  if (REMDtraj_.SetupIOarray(argIn, SetTraj().Counter(), cInfo_, Traj().Parm())) return 1;

  // Check that replica dimension valid for desired indices.
  if (targetType_ == ReplicaInfo::INDICES && 
      cInfo_.ReplicaDimensions().Ndims() != (int)remdtrajidx_.size())
  {
    mprinterr("Error: Replica # of dim (%i) not equal to target # dim (%zu)\n",
              cInfo_.ReplicaDimensions().Ndims(), remdtrajidx_.size());
    return 1;
  }

  // If target type is temperature make sure there is temperature info.
  if (targetType_ == ReplicaInfo::TEMP && !cInfo_.HasTemp()) {
    mprinterr("Error: Some or all replicas are missing temperature info.\n");
    return 1;
  }

  return 0;
}

// Trajin_Multi::BeginTraj() 
int Trajin_Multi::BeginTraj() {
  // Open the trajectories
  if (debug_ > 0) mprintf("\tREMD: OPENING %zu REMD TRAJECTORIES\n", REMDtraj_.size());
  for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
  {
    if ( (*rep)->openTrajin()) {
      mprinterr("Error: Could not open replica # %zu, '%s'\n",
                rep - REMDtraj_.begin(), REMDtraj_.f_name(rep-REMDtraj_.begin()));
      return 1;
    }
  }
  // Initialize counter.
  SetTraj().Counter().Begin();
  return 0;
}

// Trajin_Multi_EndTraj()
void Trajin_Multi::EndTraj() {
  for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
    (*rep)->closeTraj();
}

// Trajin_Multi::ReadTrajFrame()
int Trajin_Multi::ReadTrajFrame( int currentFrame, Frame& frameIn ) {
  bool replicaFound = false;

  if ( targetType_ == ReplicaInfo::TEMP ) {
    // Find target temperature, exit loop when found.
    for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
    {
      // Locate the target temp/indices out of all the replicas
      if ( (*rep)->readFrame(currentFrame, frameIn)) return 1;
      if ( frameIn.Temperature() == remdtrajtemp_ ) { replicaFound = true; break; }
    }
  } else { // Indices
    // Find target indices.
    for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
    {
      replicaFound = true;
      if ( (*rep)->readFrame(currentFrame, frameIn)) return 1;
      Frame::RemdIdxType::const_iterator tgtIdx = frameIn.RemdIndices().begin();
      for (RemdIdxType::const_iterator idx = remdtrajidx_.begin(); 
                                       idx != remdtrajidx_.end(); 
                                     ++idx, ++tgtIdx)
      {
        if ( *tgtIdx != *idx ) { replicaFound = false; break; }
      }
      if (replicaFound) break;
    }
  }

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
          Traj().Filename().base());
  if (showExtended == 1) Traj().Counter().PrintFrameInfo();
  mprintf("\n");
  if (debug_ > 0) REMDtraj_.PrintIOinfo();
  if (remdtrajidx_.empty())
    mprintf("\tLooking for frames at %.2lf K\n",remdtrajtemp_);
  else {
    mprintf("\tLooking for indices [");
    for (RemdIdxType::const_iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
      mprintf(" %i", *idx);
    mprintf(" ]\n");
  }
}
