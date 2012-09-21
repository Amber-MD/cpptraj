// RemdTraj
#include <cstring> // memcpy
#include <locale>
#include <sstream>
#include "RemdTraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
RemdTraj::RemdTraj() :
  remdtrajtemp_(0.0),
  remd_indices_(0),
  lowestRepnum_(0),
  targetType_(TEMP),
  // Used for writing replica trajectories
  hasTrajout_(false),
  remdX_(0),
  remdV_(0),
  remdT_(0.0),
  remdN_(0)
{
  remdbox_[0]=0.0;
  remdbox_[1]=0.0;
  remdbox_[2]=0.0;
  remdbox_[3]=0.0;
  remdbox_[4]=0.0;
  remdbox_[5]=0.0;
}

// DESTRUCTOR
RemdTraj::~RemdTraj() {
  std::vector<TrajectoryIO*>::iterator replica;
  //fprintf(stderr,"RemdTraj Destructor.\n");
  for (replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); replica++)
    delete *replica;
  for (replica=REMDtrajout_.begin(); replica!=REMDtrajout_.end(); replica++)
    delete *replica;
  if (remdX_!=0) delete[] remdX_;
  if (remdV_!=0) delete[] remdV_;
  if (remd_indices_!=0) delete[] remd_indices_;
}

// RemdTraj::SetTargetTemp()
void RemdTraj::SetTargetTemp(double tempIn) {
  remdtrajtemp_ = tempIn;
  targetType_ = TEMP;
}

// RemdTraj::SetTargetIdx
void RemdTraj::SetTargetIdx(RemdIdxType const& idxIn) {
  remdtrajidx_ = idxIn;
  targetType_ = INDICES;
  remd_indices_ = new int[ remdtrajidx_.size() ];
}

// RemdTraj::AddReplicaTrajin()
void RemdTraj::AddReplicaTrajin(TrajectoryIO *trajio) {
  if (trajio==NULL) {
    mprinterr("Internal Error: RemdTraj::AddReplicaTrajin: TrajectoryIO is NULL.\n");
    return;
  }
  REMDtraj_.push_back( trajio );
}

// RemdTraj::AddReplicaTrajout()
void RemdTraj::AddReplicaTrajout(TrajectoryIO *trajio) {
  if (trajio==NULL) {
    mprinterr("Internal Error: RemdTraj::AddReplicaTrajout: TrajectoryIO is NULL.\n");
    return;
  }
  REMDtrajout_.push_back( trajio );
}

// RemdTraj::LowestReplicaName()
const char *RemdTraj::LowestReplicaName() {
  return FileName_.c_str();
}

// RemdTraj::Nreplicas()
int RemdTraj::Nreplicas() {
  return (int)REMDtraj_.size();
}

// RemdTraj::SetupTemperatureList()
/** For each REMD input traj get the temperature of the first frame to
  * create a complete list of replica temperatures. This list is used
  * when writing out temperature trajectories from input replica
  * trajectories.
  */
int RemdTraj::SetupTemperatureList(int natom) {
  std::vector<TrajectoryIO*>::iterator replica;
  int err = 0;
  double tempIn;

  if (remdX_!=0) delete[] remdX_;
  if (remdV_!=0) delete[] remdV_;
  // Allocate space for reading in coords and velocities
  remdN_ = natom * 3; 
  remdX_ = new double[ remdN_ ];
  if (remdX_==0) return 1;
  if (hasVelocity_) {
    remdV_ = new double[ remdN_ ];
    if (remdV_==0) return 1;
  }

  // Get a list of all temperatures present in input REMD trajectories
  TemperatureList_.clear();
  int repnum = 0;
  for (replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); replica++) {
    err=0;
    if ( (*replica)->openTraj() ) {
      err = 1;
    } else {
      // Read 1 frame to get temperature
     
      err = ((*replica)->readFrame(0,remdX_, remdV_, remdbox_, &tempIn));
      mprintf("      Rep %i T = %6.2lf\n",repnum++,tempIn);
      TemperatureList_.push_back(tempIn);
      // Close input traj
      (*replica)->closeTraj();
    }
    // If err!=0 an error occurred, bail out.
    if (err!=0) break;
  }
  if (err==0) hasTrajout_=true;
  return err;
}

// RemdTraj::PrintNoExtError()
void RemdTraj::PrintNoExtError() {
  mprinterr("Error: Traj %s has no numerical extension, required for automatic\n",BaseFileStr());
  mprinterr("       detection of replica trajectories. Expected filename format is\n");
  mprinterr("       <Prefix>.<#> (with optional compression extension, examples:\n");
  mprinterr("       Rep.traj.nc.000,  remd.x.01.gz etc.\n");
}

// RemdTraj::SearchForReplicas()
/** Assuming this RemdTraj has been set up as the lowest replica, search for
  * all other replica names set basic assuming a nameing scheme of 
  * <PREFIX>.<EXT>[.<CEXT>], where <EXT> is a numerical extension and <CEXT>
  * is an optional compression extension. 
  * \return Found replica filenames, or an empty list on error. 
  */
std::vector<std::string> RemdTraj::SearchForReplicas() {
  std::vector<std::string> ReplicaNames;
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
    mprintf("\tREMDTRAJ: FileName=[%s]\n",FileName_.c_str());
  size_t found = FileName_.find_last_of(".");
  // NOTE: Eventually allow just numbers, no prefix? 000 002 etc
  if (found == std::string::npos) {
    PrintNoExtError();
    return ReplicaNames;
  }
  // If file is not compressed, this should be the numeric extension
  if (compressType_==NO_COMPRESSION) {
    Prefix = FileName_.substr(0, found);
    ReplicaExt = FileName_.substr(found+1);
  } else {
  // If file is compressed, this should be the compression extension
  // include the dot.
    CompressExt = FileName_.substr(found);
    // Get prefix and numerical extension, without the Compression extension
    std::string compressPrefix = FileName_.substr(0,found);
    found = compressPrefix.find_last_of(".");
    if (found == std::string::npos) {
      PrintNoExtError();
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

  // Search for and add all replicas higher than this.
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

// RemdTraj::GetTemperatureName()
/** If the temperature list has been set up, return a filename with the 
  * specified replica numbers temperature appended to it. Note that the
  * replica numbers in this case start from 0 as opposed to incoming
  * replica filenames, where the lowest replica number can be arbitrary.
  */
int RemdTraj::GetTemperatureName(std::string &tname, const char *prefix, int repnum) {
  if (TemperatureList_.empty()) {
    mprinterr("Error: RemdTraj::GetTemperatureName: TemperatureList not set.\n");
    return 1;
  }
  if (prefix==0) return 1;
  if (repnum<0 || repnum >= (int)REMDtraj_.size()) return 1;

  std::string Prefix( prefix );
  std::ostringstream Temp;
  //Temp.width(6);
  Temp.precision(2);
  Temp << std::left << std::fixed << TemperatureList_[repnum];

  tname = Prefix + "." + Temp.str();

  return 0;
}

// RemdTraj::openTraj()
/** Open each trajectory in the list. */
int RemdTraj::openTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // DEBUG
  mprintf("REMD: OPENING REMD TRAJECTORIES\n");
  // Open the trajectory
  for (replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); replica++)
    if ((*replica)->openTraj()) return 1;

  // Open output trajectories for writing
  if (hasTrajout_) {
    for (replica=REMDtrajout_.begin(); replica!=REMDtrajout_.end(); replica++) 
      if ((*replica)->openTraj()) return 1;
  }

  return 0;
}

// RemdTraj::closeTraj()
/** Close all trajectories in the list. */
void RemdTraj::closeTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // Close the trajectory
  for (replica=REMDtraj_.begin(); replica!=REMDtraj_.end(); replica++)
    (*replica)->closeTraj();

  // Close output trajectories
  if (hasTrajout_) {
    for (replica=REMDtrajout_.begin(); replica!=REMDtrajout_.end(); replica++)
      (*replica)->closeTraj();
  }
}

bool RemdTraj::IsTarget(double tempIn) {
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

// RemdTraj::readFrame()
/** Read the next frame from the list of input trajectories. Choose the
  * one that matches remdtrajtemp. Write trajectories if specified.
  */
int RemdTraj::readFrame(int set, double *X, double *V, double *box, double *T) {
  std::vector<TrajectoryIO *>::iterator reptrajin;
  std::vector<TrajectoryIO *>::iterator reptrajout;
  int trajout, nrep;
  bool found = false;
  
  for (reptrajin=REMDtraj_.begin(); reptrajin!=REMDtraj_.end(); ++reptrajin) {
    // No conversion to replica trajectories: Just find target temp
    if (!hasTrajout_) {
      if ((*reptrajin)->readFrame(set,X,V,box,T)) return 1;
      (*reptrajin)->readIndices(set,remd_indices_);
      // Check if this is the target temp
      if ( IsTarget(*T) ) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        //mprintf("REMDTRAJ: Replica %zu matches!\n", reptrajin - REMDtraj_.begin());
        return 0;
      }

    // All input REMD trajectories converted to temperature trajectories. 
    } else {
      // remdX, remdV, remdN is allocated in SetupTemperatureList 
      if ((*reptrajin)->readFrame(set,remdX_,remdV_,remdbox_,&remdT_)) return 1;
      (*reptrajin)->readIndices(set,remd_indices_);
      // Check if this is the target temp. If so, set main Frame coords/box/temp
      if ( IsTarget(remdT_) ) {
        //mprintf("REMDTRAJ: remdout: Set %i TEMP=%lf\n",set+1,remdT);
        memcpy(X, remdX_, remdN_ * sizeof(double));
        if (V!=0 && hasVelocity_) 
          memcpy(V, remdV_, remdN_*sizeof(double));
        memcpy(box, remdbox_, 6 * sizeof(double));
        *T = remdT_;
        found=true;
      }
      // Figure out which output file matches this temperature
      trajout=-1; // Will hold the index of target output traj for T
      nrep=0;
      for (reptrajout=REMDtrajout_.begin(); reptrajout!=REMDtrajout_.end(); reptrajout++) {
        if (TemperatureList_[nrep] == remdT_) {
          trajout = nrep;
          break;
        }
        ++nrep;
      }
      if (trajout==-1) {
        mprinterr("\nError: RemdTraj: remdout: Temperature %6.2lf not found in Temperature list.\n",
                  remdT_);
        return 1;
      }
      // Write to output traj
      //mprintf("REMDTRAJ: Write set %i, T %6.2lf, to %s\n",set+1, remdT, 
      //        (*reptrajout)->Filename());
      (*reptrajout)->writeFrame(set, remdX_,remdV_,remdbox_,remdT_);
    }
  }  // END LOOP over input remd trajectories
  if (found) return 0;
  // If we have made it here this means target was not found
  if (targetType_ == TEMP) {
    mprinterr("\nError: RemdTraj: Final repTemp value read= %lf, set %i\n",*T,set+OUTPUTFRAMESHIFT);
    mprinterr("Could not find target %lf in any of the replica trajectories.\n",
              remdtrajtemp_);
  } else {
    mprinterr("Error: RemdTraj: Could not find target indices ");
    for (RemdIdxType::iterator idx = remdtrajidx_.begin(); idx!=remdtrajidx_.end(); ++idx)
      mprinterr("%i,",*idx);
    mprinterr(" in any of the replica trajectories.\n");
  } 
  mprinterr("Check that all replica trajectory files were found and that\n");
  mprinterr("none of the trajectories are corrupted (e.g. missing a temperature etc).\n");
  return 1;
}

// RemdTraj::info()
void RemdTraj::info() {
  mprintf("REMD trajectories (%u total)\n", REMDtraj_.size());
  if (hasTrajout_) {
    mprintf("  remdout: trajectories will be converted to temperature trajectories (%u total):\n",
            REMDtrajout_.size());
    for (unsigned int repnum=0; repnum < REMDtraj_.size(); repnum++)
      mprintf("\t\t[%6.2lf]\n",TemperatureList_[repnum]);
  }
  if (remdtrajidx_.empty())
    mprintf("\tLooking for frames at %.2lf K",remdtrajtemp_);
  else {
    mprintf("\tLooking for indices [");
    for (RemdIdxType::iterator idx = remdtrajidx_.begin(); idx != remdtrajidx_.end(); ++idx)
      mprintf(" %i", *idx);
    mprintf(" ]");
  }
}
 
