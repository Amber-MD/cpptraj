#include <locale> // isdigit
#include <sstream>
#include "Trajin_Multi.h"
#include "CpptrajStdio.h"

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
std::vector<std::string> Trajin_Multi::SearchForReplicas(bool isCompressed) {
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
  std::vector<std::string> replica_filenames;

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

  // Check if replica trajectories are explicitly listed
  ArgList remdtraj_list( argIn->GetStringKey("trajnames"), "," );
  if (remdtraj_list.Nargs()==0) {
    // Automatically scan for additional REMD traj files.
    replica_filenames = SearchForReplicas(baseFile.IsCompressed());
  } else {
    // Get filenames from args of remdtraj_list
    replica_filenames.push_back( tnameIn );
    for (ArgList::const_iterator fname = remdtraj_list.begin();
                                 fname != remdtraj_list.end(); ++fname) 
       replica_filenames.push_back( *fname );
  }

  // Loop over all filenames in replica_filenames 
  bool lowestRep = true;
  bool repBoxInfo = false;
  for (std::vector<std::string>::iterator repfile = replica_filenames.begin();
                                          repfile != replica_filenames.end(); ++repfile)
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
      if (StartStopOffset( replica0, argIn )) return 1;
      // Check how many frames will actually be read
      if (setupFrameInfo() == 0) return 1;
      // If lowest rep has box info, all others must have box info.
      repBoxInfo = replica0->HasBox();
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
    }
    // Check for temperature information
    if ( !replica0->HasTemperature()) {
      mprinterr("Error: RemdTraj: Replica %s does not have temperature info.\n",
                (*repfile).c_str());
      return 1;
    }
    // Check that replica dimension valid for desired indices.
    if (targetType_ == INDICES && replica0->NreplicaDimensions() != (int)remdtrajidx_.size())
    {
      mprinterr("Error: RemdTraj: Replica %s # of dim (%i) not equal to target # dim (%zu)\n",
                (*repfile).c_str(), replica0->NreplicaDimensions(), remdtrajidx_.size());
      return 1;
    }
    lowestRep = false;
  }

  return 0;
}
