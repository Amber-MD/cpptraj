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
  REMDtraj_.SetDebug( debug_ );
  // Check for deprecated args
  if (argIn.hasKey("remdout")) {
    mprinterr("%s", TrajIOarray::DEPRECATED_remdout);
    return 1;
  }
  targetType_ = ReplicaInfo::NONE;
  // Process REMD-specific arguments
  std::vector<double> remdtrajval;
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
  } else if (argIn.Contains("remdtrajvalues")) {
    // Looking for target values
    ArgList valuesArg(argIn.GetStringKey("remdtrajvalues"), ",");
    if (valuesArg.empty()) {
      mprinterr("Error: remdtrajvalues expects comma-separated list of target values in each\n"
                "Error: dimension, '<dim1 val>,<dim2 val>,...,<dimN val>`.\n");
      return 1;
    }
    for (ArgList::const_iterator arg = valuesArg.begin(); // TODO validDouble?
                                 arg != valuesArg.end(); ++arg)
      remdtrajval.push_back( convertToDouble( *arg ) );
  } else if (argIn.Contains("remdtrajtemp")) {
    // Looking for target temperature
    remdtrajtemp_ = argIn.getKeyDouble("remdtrajtemp",0.0);
    targetType_ = ReplicaInfo::TEMP;
  }
  // Set up replica file names.
  if (REMDtraj_.SetupReplicaFilenames( tnameIn, argIn )) return 1;

  // Set up TrajectoryIO classes for all file names.
  if (REMDtraj_.SetupIOarray(argIn, SetTraj().Counter(), cInfo_, Traj().Parm())) return 1;

  // Additional setup for replica values.
  if (cInfo_.UseRemdValues()) {
    if (!remdtrajval.empty()) {
      // 'remdtrajvalues' specified.
      if (!remdtrajidx_.empty()) {
        mprinterr("Error: Specify either 'remdtrajvalues' or 'remdtrajidx', not both.\n");
        return 1;
      }
      if (cInfo_.HasReplicaDims()) {
        // Multiple dimensions.
        ReplicaDimArray const& Rdim = cInfo_.ReplicaDimensions();
        if (Rdim.Ndims() != (int)remdtrajval.size()) {
          mprinterr("Error: Replica # of dim (%i) not equal to target # dim (%zu)\n",
                    Rdim.Ndims(), remdtrajval.size());
          return 1;
        }
        // Get initial values from each trajectory.
        std::vector< std::vector<double> > DimValues( Rdim.Ndims() );
        Frame frameIn;
        frameIn.SetupFrameV( tparmIn->Atoms(), cInfo_ );
        if (BeginTraj()) return 1;
        for (TrajIOarray::const_iterator rep = REMDtraj_.begin(); rep != REMDtraj_.end(); ++rep)
        {
          if ( (*rep)->readFrame(0, frameIn)) return 1;
          for (int dim = 0; dim != Rdim.Ndims(); dim++) {
            switch (Rdim.DimType(dim)) {
              case ReplicaDimArray::TEMPERATURE :
                DimValues[dim].push_back( frameIn.Temperature() ); break;
              case ReplicaDimArray::PH          :
                DimValues[dim].push_back( frameIn.pH()          ); break;
              case ReplicaDimArray::REDOX       :
                DimValues[dim].push_back( frameIn.RedOx()       ); break;
              case ReplicaDimArray::RXSGLD      :
                DimValues[dim].push_back( frameIn.Temperature() ); break;
              default:
                mprinterr("Error: Dimension '%s' not supported by remdtrajvalues.\n"
                          "Error: Use 'remdtrajidx' instead.\n", Rdim.Description(dim));
                return 1;
            }
          }
        }
        EndTraj();
        // Determine index of target value in each dimension
        for (int dim = 0; dim != Rdim.Ndims(); dim++)
        {
          ReplicaInfo::Map<double> dMap;
          if (debug_ > 0) {
            mprintf("DEBUG: Initial %s values:", Rdim.Description(dim));
            for (std::vector<double>::const_iterator it = DimValues[dim].begin();
                                                     it != DimValues[dim].end(); ++it)
              mprintf(" %g", *it);
            mprintf("\n");
          }
          dMap.CreateMap(DimValues[dim], false);
          int tIdx = dMap.FindIndex( remdtrajval[dim] );
          mprintf("\tTarget index for '%s' value %g is %i\n", Rdim.Description(dim),
                  remdtrajval[dim], tIdx+1);
          if (tIdx < 0) {
            mprinterr("Error: Value %g not found in input ensemble for '%s' dimension (%i)\n",
                      remdtrajval[dim], Rdim.Description(dim), dim);
            mprinterr("Error: Potential values are:");
            for (ReplicaInfo::Map<double>::const_iterator it = dMap.begin();
                                                          it != dMap.end(); ++it)
              mprinterr(" %g", *it);
            mprinterr("\n");
            return 1;
          }
          // Replica index arguments start from 1
          remdtrajidx_.push_back( tIdx+1 );
        }
        targetType_ = ReplicaInfo::INDICES;
      } else {
        // Single dimension.
        // FIXME assuming temperature
        remdtrajtemp_ = remdtrajval.front();
        targetType_ = ReplicaInfo::TEMP;
      }
    }
  } else {
    // 'remdtrajvalues' not valid when values not present.
    if (!remdtrajval.empty()) {
      mprinterr("Error: 'remdtrajvalues' specified but trajectories do not contain replica values.\n");
      return 1;
    }
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
  }
  if (targetType_ == ReplicaInfo::NONE) {
    mprinterr("Error: With 'remdtraj' must specify either 'remdtrajtemp',\n"
              "Error: 'remdtrajidx', or 'remdtrajvalues'\n");
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
  } else {
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
    mprinterr("Error: Target replica not found, frame %i.\n", currentFrame+1);
    mprinterr("Error: Check that all replica trajectories were found and that\n");
    if (targetType_ == ReplicaInfo::TEMP)
      mprinterr("Error:  target temperature is valid");
    else
      mprinterr("Error:  target indices are valid");
    mprinterr(" for this ensemble.\n");
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
