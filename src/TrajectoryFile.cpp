// TrajectoryFile
#include <cstdio> // sscanf
#include <cstring> // strncmp, strlen
#include <cctype> // isdigit, isspace
#include <sstream> // for NumberedWrite
#include "TrajectoryFile.h"
#include "PDBfile.h" // TODO: Remove when ID_traj introduced
#include "Mol2File.h" // ditto
#include "NetcdfFile.h" // ditto
#include "CpptrajStdio.h"
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
#ifdef BINTRAJ
  #include "Traj_AmberNetcdf.h"
  #include "Traj_AmberRestartNC.h"
#endif
#include "Traj_PDBfile.h"
#include "Traj_AmberRestart.h"
#include "Traj_Mol2File.h"
#include "Traj_Conflib.h"
#include "Traj_CharmmDcd.h"
#include "RemdTraj.h"

// CONSTRUCTOR
TrajectoryFile::TrajectoryFile() {
  debug_ = 0;
  progress_ = NULL;
  trajio_ = NULL;
  trajName_ = NULL;
  trajParm_ = NULL;
  fileAccess_ = READTRAJ;
  numFramesProcessed_ = 0;
  start_ = 0;
  stop_ = -1;
  offset_ = 1;
  total_frames_ = 0;
  total_read_frames_ = -1;
  currentFrame_ = 0;
  targetSet_ = 0;
  frameskip_ = 1;
  FrameRange_ = NULL;
  nobox_ = false;
  // NOTE: Using this instead of trajio->IsOpen since trajio doesnt always
  //       use the FileIO routines (e.g. netcdf files).
  trajIsOpen_ = false;
}

// DESTRUCTOR
TrajectoryFile::~TrajectoryFile() {
  if (trajio_!=NULL) {
    if (trajIsOpen_) EndTraj();
    delete trajio_;
  }
  if (FrameRange_!=NULL) delete FrameRange_;
  if (progress_!=NULL) delete progress_;
}

// TrajectoryFile::SetDebug()
/** Set debug level. */
void TrajectoryFile::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("  TrajectoryFile debug level set to %i\n",debug_);
}

// TrajectoryFile::ID_TrajFormat()
TrajectoryFile::TrajFormatType TrajectoryFile::ID_TrajFormat(TrajectoryIO &trajectory) 
{
  unsigned char magic[3];
  const int BUF_SIZE = 83;
  char buffer1[BUF_SIZE];
  char buffer2[BUF_SIZE];
  // TODO: Remove the next 3, ID should be in individual TrajectoryIO classes
  PDBfile pdbfile;
  Mol2File mol2file;
  NetcdfFile ncfile;

  // Open trajectory
  if (trajectory.IsOpen())
    trajectory.CloseFile();
  if (trajectory.OpenFile())
    return UNKNOWN_TRAJ;
  // Get magic number
  magic[0] = 0; magic[1] = 0; magic[2] = 0;
  // NOTE: Replace with 1 read
  trajectory.IO->Read(magic  ,1);
  trajectory.IO->Read(magic+1,1);
  trajectory.IO->Read(magic+2,1);
  trajectory.CloseFile();
  if (debug_>0) 
    mprintf("\tTrajectory Hex sig: %x %x %x", magic[0],magic[1],magic[2]);

  // NETCDF
  if (magic[0]==0x43 && magic[1]==0x44 && magic[2]==0x46) {
    if (debug_>0) mprintf("  NETCDF file\n");
    if (trajectory.IsCompressed()) {
      mprinterr("Error: Compressed NETCDF files are not currently supported.\n");
      return UNKNOWN_TRAJ;
    }
    // Determine whether this is a trajectory or restart Netcdf from the Conventions
    NetcdfFile::NCTYPE nctype = ncfile.GetNetcdfConventions( trajectory.Name() );
    switch (nctype) {
      case NetcdfFile::NC_AMBERTRAJ   : return AMBERNETCDF; break;
      case NetcdfFile::NC_AMBERRESTART: return AMBERRESTARTNC; break;
      case NetcdfFile::NC_UNKNOWN:
      default: return UNKNOWN_TRAJ;
    }
    // Shouldnt get here, just for safetys sake
    return UNKNOWN_TRAJ;
  }

  // ---------- ID by file characteristics; read the first two lines ----------
  // Initialize buffers to NULL
  memset(buffer1,' ',BUF_SIZE);
  memset(buffer2,' ',BUF_SIZE);
  buffer1[0]='\0';
  buffer2[0]='\0';
  if (trajectory.OpenFile())
    return UNKNOWN_TRAJ;
  trajectory.IO->Gets(buffer1, BUF_SIZE);
  trajectory.IO->Gets(buffer2, BUF_SIZE);
  trajectory.CloseFile();
  
  // If both lines have PDB keywords, assume PDB
  trajectory.OpenFile();
  if ( pdbfile.ID( trajectory.IO ) )
  {
    if (debug_>0) mprintf("  PDB file\n");
    trajectory.CloseFile();
    return PDBFILE;
  }

  // If either buffer contains a TRIPOS keyword assume Mol2
  // NOTE: This will fail on tripos files with extensive header comments.
  //       A more expensive check for mol2 files is below.
  trajectory.OpenFile();
  if ( mol2file.ID( trajectory.IO ) )
  {
    if (debug_>0) mprintf("  TRIPOS MOL2 file\n");
    trajectory.CloseFile();
    return MOL2FILE;
  }
  trajectory.CloseFile();

  // If the second 5 chars are C O R D, assume charmm DCD
  if (strncmp(buffer1+4,"CORD",4)==0) {
    if (debug_>0) mprintf("  CHARMM DCD file\n");
    return CHARMMDCD;
  }

  // Amber Restart
  // Check for an integer (I5) followed by 0-2 scientific floats (E15.7)
  if (strlen(buffer2)<=36) {
    //mprintf("DEBUG: Checking restart.\n");
    //mprintf("DEBUG: buffer2=[%s]\n",buffer2);
    int i=0;
    for (; i<5; i++) {
      if (!isspace(buffer2[i]) && !isdigit(buffer2[i])) break;
      //mprintf("DEBUG:    %c is a digit/space.\n",buffer2[i]);
    }
    //mprintf("DEBUG: i=%i\n");
    //if ( i==5 && strchr(buffer2,'E')!=NULL ) {
    if ( i==5 ) {
      if (debug_>0) mprintf("  AMBER RESTART file\n");
      return AMBERRESTART;
    }
  }

  // Check if second line contains REMD/HREMD, Amber Traj with REMD header
  if ( strncmp(buffer2,"REMD",4)==0 ||
       strncmp(buffer2,"HREMD",5)==0   ) 
  {
    if (debug_>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
    return AMBERTRAJ;
  }

  // Check if we can read at least 3 coords of width 8, Amber trajectory
  float TrajCoord[3];
  if ( sscanf(buffer2, "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 ) 
  {
    if (debug_>0) mprintf("  AMBER TRAJECTORY file\n");
    return AMBERTRAJ;
  }
  
/*  // ---------- MORE EXPENSIVE CHECKS ----------
  // Reopen and scan for Tripos mol2 molecule section
  // 0 indicates section found.
  if (trajectory.OpenFile())
    return UNKNOWN_TRAJ;
  if (Mol2ScanTo(trajectory.IO, MOLECULE)==0) {
    if (debug_>0) mprintf("  TRIPOS MOL2 file\n");
    trajectory.CloseFile();
    return MOL2FILE;
  }
  trajectory.CloseFile();*/

  // ---------- EXPERIMENTAL ----------
  // If the file format is still undetermined and the file name is conflib.dat,
  // assume this is a conflib.dat file from LMOD. Cant think of a better way to
  // detect this since there is no magic number but the file is binary.
  if ( trajectory.FilenameIs((char*)"conflib.dat") ) 
  {
    mprintf("  LMOD CONFLIB file\n");
    return CONFLIB;
  }

  // Unidentified file
  mprintf("  Warning: %s: Unknown trajectory format.\n",trajName_);
  return UNKNOWN_TRAJ;
}

// TrajectoryFile::setupRemdTrajIO()
/** Assuming that trajio has already been set up as the lowest replica (and 
  * checked for temperature information), set this trajectory up as an REMD 
  * trajectory. A special trajio object will be set up that itself contains 
  * a list of trajio objects, each one corresponding to a replica. If a
  * filename is given via <remdout>, the replica trajectories will be converted
  * to temperature trajectories with names <remdout>.<Temperature> and format
  * remdfmt. Returns a RemdTraj trajio class.
  */
TrajectoryIO *TrajectoryFile::setupRemdTrajIO(double remdtrajtemp, char *remdout, 
                                              TrajFormatType remdfmt,
                                              ArgList &remdtraj_list) 
{
  RemdTraj *remdio=NULL;
  TrajectoryIO *replica0=NULL;
  std::string repFilename;
  std::vector<std::string> replica_filenames;
  char *fname;

  // Initial set up of remd trajio object
  remdio = new RemdTraj();
  // Make remdio a copy of trajio, which should be the lowest replica.
  remdio->TrajectoryIO::operator=( *trajio_ );
  // Set target temperature
  remdio->SetTargetTemp(remdtrajtemp);
  // Add trajio as first replica.
  remdio->AddReplicaTrajin( trajio_ ); 

  // ------------------------------------------------------
  // Automatically scan for additional REMD traj files.
  if (remdtraj_list.Nargs()==0) {
    replica_filenames = remdio->SearchForReplicas();
  // Get filenames from args of remdtraj_list
  } else {
    while ( (fname = remdtraj_list.getNextString()) != NULL ) {
       repFilename.assign( fname );
       replica_filenames.push_back( repFilename );
    }
  }

  // Loop over all filenames in replica_filenames
  int repnum = 1;
  for (std::vector<std::string>::iterator repfile = replica_filenames.begin();
                                          repfile != replica_filenames.end();
                                          repfile++)
  {
    mprintf("\t[%s]\n",(*repfile).c_str());
    // Set up replica file repnum trajectory IO
    replica0 = setupTrajIO((char*)(*repfile).c_str(),READTRAJ,UNKNOWN_TRAJ);
    if (replica0==NULL) {
      mprinterr("    Error: RemdTraj: Could not set up replica %i file %s\n",
                repnum,(*repfile).c_str());
      delete remdio;
      return NULL;
    }
    // Pushing replica0 here allows the remdio destructor to handle it on errors
    remdio->AddReplicaTrajin( replica0 );
    // Check that number of frames matches
    int repframes = replica0->setupTrajin(trajParm_);
    if (repframes < 0 || repframes != total_frames_) {
      mprinterr("    Error: RemdTraj: Replica %i frames (%i) does not match\n",
                repnum,repframes);
      mprinterr("           # frames in replica 0 (%i).\n",total_frames_);
      delete remdio;
      return NULL;
    }
    // Check for temperature information
    if ( !replica0->HasTemperature()) {
      mprinterr("    Error: RemdTraj: Replica %i does not have temperature info.\n",repnum);
      delete remdio;
      return NULL;
    }
    // Check box information
    if ( replica0->HasBox() != remdio->HasBox() ) {
      mprinterr("    Error: RemdTraj: Replica %i box info does not match lowest replica.\n");
      delete remdio;
      return NULL;
    }
    // Increment
    repnum++;
  }

  // If remdout was specified, set up output trajectories
  if (remdout!=NULL) {
    // Set up temperature list
    if ( remdio->SetupTemperatureList(trajParm_->Natom()) ) {
      mprinterr("Error: RemdTraj: remdout: could not get temperature list.\n");
      delete remdio;
      return NULL;
    }
    // BEGIN LOOP over Temperature
    for (repnum=0; repnum < remdio->Nreplicas(); repnum++) { 
      // Set up output filename for this temperature
      if (remdio->GetTemperatureName(repFilename,remdout,repnum)!=0) {
        mprinterr("Error: Could not get temperature filename for replica %i\n",repnum);
        delete remdio;
        return NULL;
      }
      mprintf("    Creating remd output traj: %s\n",repFilename.c_str());
      // Set up file with given type and format, and set up replica0 for the format.
      if ( (replica0=setupTrajIO((char*)repFilename.c_str(),WRITETRAJ,remdfmt))==NULL ) 
      {
        mprinterr("    Error: Could not set up T-replica file %s for writing.\n",
                  repFilename.c_str());
        delete remdio;
        return NULL;
      }
      // Add to the remd trajout list
      remdio->AddReplicaTrajout( replica0 );
      // Set up write here, trajParm will not change
      if (remdio->HasBox()) replica0->SetBox();
      if (replica0->setupTrajout(trajParm_)) {
        delete remdio;
        return NULL;
      }
    } // END LOOP over input remd trajectories
  } // END REMDtrajout

  // Since the repeated calls to setupTrajIO have overwritten trajName, 
  // reset it to the lowest replica name. 
  //SetTrajName( remdio->GetLowestReplicaName() );
  //SetTrajName( (char*)remdio->LowestReplicaName() );
  trajName_ = remdio->BaseName();

  return (TrajectoryIO*) remdio;
}

// TrajectoryFile::setupTrajIO()
/** Set up basic trajectory file for given access type. Set the trajectory
  * name to be the base filename. Return the trajectory IO object for 
  * the format.
  */
TrajectoryIO *TrajectoryFile::setupTrajIO(char *tname, TrajAccessType accIn, 
                                          TrajFormatType fmtIn) 
//, FileType typeIn) 
{
  TrajectoryIO *tio = NULL;
  TrajectoryIO basicTraj;
  int err = 1;
  TrajFormatType trajFormat = fmtIn;

  switch (accIn) {
    case READTRAJ  : err = basicTraj.SetupRead(tname,debug_); break;
    //case WRITE : err = basicTraj.SetupWrite(tname,fmtIn,typeIn,debug_); break;
    case WRITETRAJ : err = basicTraj.SetupWrite(tname,debug_); break;
    case APPENDTRAJ: err = basicTraj.SetupAppend(tname,debug_); break;
  }
  if (err!=0) {
    //mprinterr("    Error: Could not set up file %s.\n",tname);
    return NULL;
  }
  // Set trajectory name to the base filename
  // NOTE: Case back to char* ok?
  //SetTrajName( (char*)basicTraj.BaseName() );

  // If not writing, determine the file format
  if (accIn != WRITETRAJ)
    trajFormat = ID_TrajFormat( basicTraj );

  // Allocate trajectory IO type based on format
  switch ( trajFormat ) {
    case AMBERRESTART: tio = new AmberRestart(); break;
    case AMBERTRAJ   : tio = new AmberCoord();    break;
    case AMBERNETCDF :
#ifdef BINTRAJ
      tio = new AmberNetcdf();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      tio = new AmberRestartNC();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#endif
      break;
    case PDBFILE     : tio = new Traj_PDBfile();      break;
    case CONFLIB     : tio = new Conflib();      break;
    case MOL2FILE    : tio = new Traj_Mol2File();     break;
    case CHARMMDCD   : tio = new CharmmDcd();    break;
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",tname);
      return NULL;
  }

  // Happens when memory cannot be allocd, or not compiled for netcdf
  if (tio==NULL) return NULL;
 
  // Place the basic file in the trajectory IO class
  tio->TrajectoryIO::operator=( basicTraj );

   // Set debug level
  tio->SetDebug(debug_);

  // Set filename
  trajName_ = tio->BaseName();

  return tio;
}

// TrajectoryFile::SetArgs()
/** For reading trajectories only, called after initial trajectory setup. 
  * Set the start, stop, and offset args based on user input. Do some bounds 
  * checking.
  * For compatibility with ptraj, frame args start at 1. Internal frame #s 
  * start at 0. So for a traj with 10 frames:
  * - Internal #: 0 1 2 3 4 5 6 7 8 9
  * - Frame Arg#: 1 2 3 4 5 6 7 8 9 10
  * - Defaults: startArg=1, stopArg=-1, offsetArg=1
  */
int TrajectoryFile::SetArgs(ArgList *argIn) {
  int startArg;
  int stopArg;
  int offsetArg;

  if (argIn==NULL) return 0;
  // lastframe is a special case where only the last frame will be selected
  if (argIn->hasKey("lastframe")) {
    if (total_frames_>0) {
      startArg = total_frames_;
      stopArg = total_frames_;
      offsetArg = 1;
    } else {
      mprinterr("Error: [%s] lastframe specified but # frames could not be determined.\n",
                trajName_);
      return 1;
    }
  } else {
    startArg = argIn->getNextInteger(1);
    if (argIn->hasKey("last"))
      stopArg = -1;
    else  
      stopArg = argIn->getNextInteger(-1);
    offsetArg = argIn->getNextInteger(1);
  }

#ifdef DEBUGTRAJ
  mprintf("DEBUG [%s] SetArgs: Original start, stop: %i %i\n",trajName_,start_,stop_);
  mprintf("DEBUG [%s] SetArgs: Original startArg, stopArg: %i %i\n",trajName_,startArg,stopArg);
#endif
  if (startArg!=1) {
    if (startArg<1) {
      mprintf("    Warning: %s start argument %i < 1, setting to 1.\n",trajName_,startArg);
      start_ = 0; // cpptraj = ptraj - 1
    } else if (total_frames_>=0 && startArg>total_frames_) {
      // If startArg==stopArg and is greater than # frames, assume we want
      // the last frame (useful when reading for reference structure).
      if (startArg==stopArg) {
        mprintf("    Warning: %s start %i > #Frames (%i), setting to last frame.\n",
                trajName_,startArg,total_frames_);
        start_ = total_frames_ - 1;
      } else {
        mprinterr("Error: [%s] start %i > #Frames (%i), no frames will be processed.\n",
                trajName_,startArg,total_frames_);
        //start=startArg - 1;
        return 1;
      }
    } else
      start_ = startArg - 1;
  }
  if (stopArg!=-1) {
    if ((stopArg - 1)<start_) { // cpptraj = ptraj - 1
      mprinterr("Error: [%s] stop %i < start, no frames will be processed.\n",
              trajName_,stopArg);
      //stop = start;
      return 1;
    } else if (total_frames_>=0 && stopArg>total_frames_) {
      mprintf("    Warning: %s stop %i >= #Frames (%i), setting to max.\n",
              trajName_,stopArg,total_frames_);
      stop_ = total_frames_;
    } else
      stop_ = stopArg;
  }
  if (offsetArg!=1) {
    if (offsetArg<1) {
      mprintf("    Warning: %s offset %i < 1, setting to 1.\n",
              trajName_,offsetArg);
      offset_ = 1;
    } else if (stop_!=-1 && offsetArg > stop_ - start_) {
      mprintf("    Warning: %s offset %i is so large that only 1 set will be processed.\n",
              trajName_,offsetArg);
      offset_ = offsetArg;
    } else
      offset_ = offsetArg;
  }
  if (debug_>0)
    mprintf("DEBUG [%s] SetArgs: Start %i Stop %i  Offset %i\n",trajName_,start_,stop_,offset_);
  return 0;
}

// TrajectoryFile::SingleFrame()
/** Tell the trajectory to set up stop and offset so that only start frame
  * will be processed.
  */
void TrajectoryFile::SingleFrame() {
  stop_ = start_ + 1;
  offset_ = 1;
  // Call setupFrameInfo to recalc total_read_frames. Since setupFrameInfo 
  // should have already been called in SetupRead (and thus any errors 
  // handled there) dont check for an error here. It should return 1.
  if ( setupFrameInfo() != 1 ) {
    mprintf("  Warning: Single frame requested for %s but not calcd!\n",trajName_);
    mprintf("           start/stop/offset (%i, %i, %i)\n",start_+1,stop_+1,offset_);
  }
}

// TrajectoryFile::SetupRead()
/** Set up trajectory for reading. Input trajectory filename can be specified
  * explicitly, or if not it should be the second argument in the given
  * argument list. Associate this trajectory with the given parm file.
  */
int TrajectoryFile::SetupRead(char *tnameIn, ArgList *argIn, Topology *tparmIn) {
  char *tname = NULL;
  // REMD
  double remdtrajtemp = 0.0;
  char *remdout = NULL;
  bool remdtraj=false;
  TrajFormatType remdfmt = AMBERTRAJ;
  ArgList remdtraj_list;

  // Mark as not yet open
  trajIsOpen_ = false;

  // Check for trajectory filename
  if (tnameIn==NULL && argIn!=NULL)
    tname = argIn->getNextString();
  else
    tname = tnameIn;
  if (tname==NULL) {
    mprinterr("Error: TrajectoryFile::SetupRead: Filename is NULL.\n");
    return 1;
  }

  // Check that file exists
  if (!fileExists(tname)) {
    mprinterr("Error: File %s does not exist.\n",tname);
    return 1;
  }

  // Check for associated parm file
  if (tparmIn==NULL) {
    mprinterr("Error: TrajectoryFile::SetupRead: Parm file is NULL.\n");
    return 1;
  }
  trajParm_ = tparmIn;

  // Check for remdtraj keyword; if present, get args pertaining to REMD traj.
  if (argIn!=NULL && argIn->hasKey("remdtraj")) {
    remdtraj=true;
    // Get target temperature
    remdtrajtemp=argIn->getKeyDouble("remdtrajtemp",0.0);
    // If remdout specified, get a filename and optionally a format keyword
    // (default AMBERTRAJ if none specified) that will be used for writing 
    // out temperature trajectories. Currently only AmberNetcdf and AmberTraj
    // formats supported for this option.
    remdout = argIn->getKeyString("remdout",NULL);
    if (remdout!=NULL) 
      remdfmt = GetFormatFromArg(argIn);
    if (remdfmt!=AMBERTRAJ && remdfmt!=AMBERNETCDF) {
      mprinterr("Error: remdout (%s): Unsupported format. Currently only amber\n",remdout);
      mprinterr("       trajectory and amber netcdf files supported for remdout.\n");
      return 1;
    }
    // Check if replica trajectories are explicitly listed
    remdtraj_list = argIn->getKeyArgList("trajnames");
  }

  // Set this trajectory access to READ
  fileAccess_ = READTRAJ;

  // Set up file; among other things this will determine the type and 
  // format, and set up trajio for the format.
  if ( (trajio_ = setupTrajIO(tname,READTRAJ,UNKNOWN_TRAJ))==NULL ) {
    mprinterr("    Error: Could not set up file %s for reading.\n",tname);
    return 1;
  }

  // Set up the format for reading and get the number of frames.
  // -1 indicates an error.
  // -2 indicates the number of frames could not be determined, read to EOF.
  total_frames_ = trajio_->setupTrajin(trajParm_);
  if (total_frames_ == -1) {
    mprinterr("    Error: Could not set up %s for reading.\n",tname);
    return 1;
  }
  if (total_frames_>-1)
    mprintf("\t[%s] contains %i frames.\n",trajName_,total_frames_);
  else
    mprintf("\t[%s] contains an unknown number of frames.\n",trajName_);

  // Set stop based on calcd number of frames.
  if (total_frames_==0) {
    mprinterr("  Error: trajectory %s contains no frames.\n",trajName_);
    return 1;
  }
  if (total_frames_>0)
    stop_ = total_frames_; 
  else
    stop_ = -1;

  // If the trajectory has box coords, set the box type from the box Angles.
  if (trajio_->CheckBoxInfo(trajParm_)!=0) {
    mprinterr("Error in trajectory %s box information.\n",trajName_);
    return 1;
  }

  // Set start, stop, and offset args from user input.
  if (SetArgs(argIn)) return 1;

  // Call setupFrameInfo to calc actual start and stop values based on
  // offset, as well as total_read_frames
  if ( setupFrameInfo() == 0 ) {
    mprinterr("  Error: No frames will be read from %s based on start, stop,\n",trajName_);
    mprinterr("         and offset values (%i, %i, %i)\n",start_+1,stop_+1,offset_);
    return 1;
  }

  // For replica trajectories, replace the current trajio (which should be
  // the lowest replica) with a special trajio object that will have 
  // a list of trajio objects, each one pertaining to a different replica.
  if (remdtraj) {
    // Check that this trajio has temperature info
    if (!trajio_->HasTemperature()) {
      mprinterr("Error: RemdTraj: Lowest replica file %s does not have temperature info.\n",
                trajName_);
      return 1;
    }
    // Replace this trajio with a special replica traj one that will contain
    // trajio objects for all replicas
    if ( (trajio_ = setupRemdTrajIO(remdtrajtemp, remdout, remdfmt, remdtraj_list))==NULL )
      return 1;
  }
  
  return 0;
}

// TrajectoryFile::setupFrameInfo()
/** Calculate number of frames that will be read based on start, stop, and
  * offset (total_read_frames). 
  * \return the total number of frames that will be read for this traj.
  * \return -1 if the number of frames could not be determined.
  */
int TrajectoryFile::setupFrameInfo() {
  int Nframes;
  int ptraj_start_frame, ptraj_end_frame;
  int traj_start_frame, traj_end_frame;
  // DEBUG - No mpi for now
  int worldrank = 0;
  int worldsize = 1;

  //mprintf("DEBUG: Calling setupFrameInfo for %s with %i %i %i\n",trajName,
  //        start,stop,offset);
  if (stop_==-1) return -1;

  // Calc total frames that will be read
  Nframes = stop_ - start_;
  total_read_frames_ = Nframes / offset_;
  // Round up
  if ( (Nframes % offset_) > 0 )
    ++total_read_frames_;
  //divresult = div( (stop - start), offset);
  //total_read_frames = divresult.quot;
  //if (divresult.rem!=0) total_read_frames++;

  // Calc min num frames read by any given thread
  // last thread gets leftovers
  // In case of 0, last thread gets the frame
  Nframes = total_read_frames_ / worldsize;
  int leftover_frames = total_read_frames_ % worldsize;
  //divresult = div(total_read_frames,worldsize);
  //Nframes=divresult.quot;

  // Ptraj (local) start and end frame
  ptraj_start_frame = (worldrank*Nframes);
  ptraj_end_frame = ptraj_start_frame + Nframes;
  // Last thread gets the leftovers
  if (worldrank == worldsize-1) 
    ptraj_end_frame += leftover_frames;

  // Actual Traj start and end frame (for seeking)
  traj_start_frame = (ptraj_start_frame * offset_) + start_;
  traj_end_frame = ((ptraj_end_frame-1) * offset_) + start_;

  start_ = traj_start_frame;
  stop_ = traj_end_frame;

#ifdef TRAJDEBUG
  mprintf("DEBUG SETUPFRAMEINFO: %i-%i total %i\n",start_,stop_,total_read_frames_);
#endif

  return total_read_frames_;
}

// TrajectoryFile::GetFormatFromArg()
/** Given an arglist, search for one of the file format keywords.
  * Default to AmberTraj if no arglist given or no keywords present. 
  */
TrajectoryFile::TrajFormatType TrajectoryFile::GetFormatFromArg(ArgList *argIn) 
{
  TrajFormatType writeFormat = AMBERTRAJ;
  if (argIn==NULL) return writeFormat;
  if      ( argIn->hasKey("pdb")      ) writeFormat=PDBFILE;
  else if ( argIn->hasKey("netcdf")   ) writeFormat=AMBERNETCDF;
  else if ( argIn->hasKey("restart")  ) writeFormat=AMBERRESTART;
  else if ( argIn->hasKey("ncrestart")) writeFormat=AMBERRESTARTNC;
  else if ( argIn->hasKey("restartnc")) writeFormat=AMBERRESTARTNC;
  else if ( argIn->hasKey("mol2")     ) writeFormat=MOL2FILE;
  else if ( argIn->hasKey("dcd")      ) writeFormat=CHARMMDCD;
  return writeFormat;
}

std::string TrajectoryFile::GetExtensionForType(TrajFormatType typeIn) {
  std::string ext;
  switch (typeIn) {
    case PDBFILE : ext=".pdb"; break;
    case AMBERTRAJ: ext=".crd"; break;
    case AMBERNETCDF: ext=".nc"; break;
    case AMBERRESTART: ext=".rst7"; break;
    case AMBERRESTARTNC: ext=".ncrst"; break;
    case MOL2FILE: ext=".mol2"; break;
    case CHARMMDCD: ext=".dcd"; break;
    default: ext="";
  }
  return ext;
}
  

// TrajectoryFile::SetupWriteWithArgs()
/** Like SetupWrite, but intended for internal use. Allows a static
  * space-separated string to be passed in, which will be converted
  * to an argument list and passed to SetupWrite.
  */
int TrajectoryFile::SetupWriteWithArgs(char *tnameIn, const char *argstring,
                                       Topology *tparmIn, TrajFormatType fmtIn) 
{
  ArgList tempArg;
  //char *tempString;
  //tempArg.SetDebug(2);
  // Since ArgList uses strtok cannot pass it a const string, copy to 
  // temporary string.
  //tempString = new char[ strlen(argstring) + 1];
  //strcpy(tempString,argstring);
  tempArg.SetList((char*)argstring, " ");
  //delete[] tempString;
  return SetupWrite(tnameIn,&tempArg,tparmIn,fmtIn);
}

// TrajectoryFile::SetupWrite()
// Currently only used byt Clustering
int TrajectoryFile::SetupWrite(char *tnameIn, Topology *tparmIn, char *fmtArg) {
  ArgList tempArg;
  if (fmtArg!=NULL)
    tempArg.AddArg(fmtArg);
  return SetupWrite(tnameIn, &tempArg, tparmIn, UNKNOWN_TRAJ);
}

// TrajectoryFile::SetupNumberedWrite()
// CUrrently only used by clustering
int TrajectoryFile::SetupNumberedWrite(char *prefixIn, int numIn, Topology *tparmIn, 
                                       char *fmtArg)
{
  std::string Prefix;
  std::string Ext;
  std::ostringstream Num;
  ArgList trajFormatArg;

  if (prefixIn==NULL) return 1;
  Prefix.assign(prefixIn);
  if (!(Num << numIn)) return 1;
  Prefix += ".";
  Prefix += Num.str();
  if (fmtArg!=NULL)
    trajFormatArg.AddArg(fmtArg);
  TrajFormatType trajFmt = GetFormatFromArg(&trajFormatArg);
  Ext = GetExtensionForType(trajFmt);
  Prefix += Ext;
  return SetupWrite((char*)Prefix.c_str(), NULL, tparmIn, trajFmt);
}

// TrajectoryFile::SetupWrite()
/** Set up trajectory for writing. Output trajectory filename can be specified
  * explicitly, or if not it should be the second argument in the given
  * argument list. Associate with the given parm file initially, but setup for 
  * the format is performed on the first write call to accomodate state changes
  * like stripped atoms and so on. 
  * NOTE: make remdtraj a generic trigger for hasTemperature?
  */
int TrajectoryFile::SetupWrite(char *tnameIn, ArgList *argIn, Topology *tparmIn,
                               TrajFormatType writeFormatIn) 
{
  char *tname = NULL;
  fileAccess_ = WRITETRAJ;
  TrajFormatType writeFormat = writeFormatIn;;
  char *onlyframes=NULL;

  // Mark as not yet open 
  trajIsOpen_ = false;

  // Check for trajectory filename
  if (tnameIn==NULL && argIn!=NULL)
    tname = argIn->getNextString();
  else
    tname = tnameIn;
  if (tname==NULL) {
    mprinterr("Error: TrajectoryFile::SetupWrite: Filename is NULL.\n");
    return 1;
  }

  // Check for associated parm file
  if (tparmIn==NULL) {
    mprinterr("Error: TrajectoryFile::SetupWrite: Parm file is NULL.\n");
    return 1;
  }
  trajParm_ = tparmIn;

  // If a write format was not specified (UNKNOWN_TRAJ) check the argument
  // list to see if format was specified there. Defaults to AMBERTRAJ.
  if (writeFormat==UNKNOWN_TRAJ)
    writeFormat = GetFormatFromArg(argIn);

  // Check for append keyword
  if (argIn!=NULL && argIn->hasKey("append")) fileAccess_ = APPENDTRAJ;

  // Set up file with given type and format, and set up trajio for the format.
  if( (trajio_ = setupTrajIO(tname,fileAccess_,writeFormat))==NULL ) {
    mprinterr("    Error: Could not set up file %s for writing.\n",tname);
    return 1;
  }

  // Process additional arguments
  if (argIn!=NULL) {
    // Get specified title if any - will not set if NULL
    trajio_->SetTitle( argIn->getKeyString("title", NULL) );

    // Get a frame range for trajout
    onlyframes = argIn->getKeyString("onlyframes",NULL);
    if (onlyframes!=NULL) {
      FrameRange_ = new Range();
      if ( FrameRange_->SetRange(onlyframes) ) {
        mprintf("Warning: trajout %s: onlyframes: %s is not a valid range.\n",tname,onlyframes);
        delete FrameRange_;
      } else {
        FrameRange_->PrintRange("      Saving frames",0);
      }
      // User frame args start from 1. Start from 0 internally.
      FrameRange_->ShiftBy(-1);
    }

    // Check for nobox argument - will override any box info present in parm
    // when trajectory IO is set up.
    nobox_ = argIn->hasKey("nobox");

    // Process any write arguments specific to certain formats not related
    // to parm file. Options related to parm file are handled on the first
    // write in WriteFrame.
    if (trajio_->processWriteArgs(argIn)) {
      mprinterr("Error: trajout %s: Could not process arguments.\n",tname);
      return 1;
    }
  }

  // No more setup here; Write is set up when first frame written.
  return 0;
}

// TrajectoryFile::BeginTraj()
/** Prepare trajectory file for reading/writing. For reads, open the traj and 
  * set up a progress bar if requested. Not needed for writes since opening
  * occurs when WriteFrame() is called for the first time.
  */
int TrajectoryFile::BeginTraj(bool showProgress) {
  // For writes the trajectory is opened on first write in WriteFrame 
  if (fileAccess_ != READTRAJ) return 0;

  // Open the trajectory
  if (trajio_->openTraj()) {
    mprinterr("Error: TrajectoryFile::BeginTraj: Could not open %s\n",trajName_);
    return 1;
  }
  trajIsOpen_ = true;
  numFramesProcessed_ = 0;

  // Set up a progress bar
  if (showProgress) progress_ = new ProgressBar(total_read_frames_);

  // Determine what frames will be read
  targetSet_ = start_;
  if (trajio_->Seekable()) {
    frameskip_ = offset_;
    currentFrame_ = start_;
  } else {
    frameskip_ = 1;
    currentFrame_ = 0;
  }

  return 0;
}

// TrajectoryFile::PrintInfoLine()
/** Print a smaller amount of information than the PrintInfo function to
  * a single line. Give the trajectory name, start, stop, and offset.
  */
void TrajectoryFile::PrintInfoLine() {
  //rprintf( "----- [%s] (%i-%i, %i) -----\n",trajName,currentFrame+1,stop+1,offset);
  if (stop_ != -1)
    rprintf( "----- [%s] (%i-%i, %i) -----\n",trajName_,start_+1,stop_+1,offset_);
  else
    rprintf( "----- [%s] (%i-EOF, %i) -----\n",trajName_,start_+1,offset_);
}

// TrajectoryFile::EndTraj()
/** Close the trajectory. Should be called after all reads / writes completed.
  */
int TrajectoryFile::EndTraj() {
  trajio_->closeTraj();
  trajIsOpen_ = false;
  return 0;
}

// TrajectoryFile::GetNextFrame()
/** Get the next target frame from trajectory. Update the number of frames
  * read while getting to target (if traj is seekable this will always be 1).
  * \return 1 on successful read.
  * \return 0 if no frames could be read.
  */
// NOTE: Just check if traj is finished after targetSet is incremented?
int TrajectoryFile::GetNextFrame(Frame& FrameIn) { 
  bool tgtFrameFound;
#ifdef TRAJDEBUG
  mprinterr("Getting frame %i from %s (stop=%i)\n",currentFrame_,trajName_,stop_);
#endif
  // If the current frame is out of range, exit
  if (currentFrame_>stop_ && stop_!=-1) return 0;
  //if (currentFrame>stop) return 0;
  if (progress_!=NULL) 
    progress_->Update(numFramesProcessed_);
    //progress->PrintBar(currentFrame);

  tgtFrameFound=false;

  while ( !tgtFrameFound ) {
#ifdef TRAJDEBUG
    mprinterr("Attempting read of frame %i from %s\n",currentFrame_,trajName_);
#endif
    if (trajio_->readFrame(currentFrame_,
                           FrameIn.X_,FrameIn.V_,FrameIn.box_,&(FrameIn.T_))) 
      return 0;
    //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
#ifdef TRAJDEBUG
    mprinterr("Frame %i has been read from %s (target=%i)\n",currentFrame_,trajName_,targetSet_);
#endif
    if (currentFrame_ == targetSet_) {
      tgtFrameFound=true;
      targetSet_ += offset_;
    }
    ++numFramesProcessed_;
    currentFrame_ += frameskip_;
  }

  return 1;
}

// TrajectoryFile::WriteFrame()
/** Write the given coordinates, box, and Temperature. Only write if the given
  * parm matches the parm that SetupWrite was called with. This allows for
  * things like stripped topologies (a stripped top will have the same pindex
  * as the original top). 
  * The first time this is called with the correct parm the output trajectory 
  * will be set up based on the parm it is called with and the file will be 
  * opened, so there is no need to call BeginTraj for output trajectories.
  * EndTraj() should still be called.
  */
int TrajectoryFile::WriteFrame(int set, Topology *tparmIn, Frame &FrameOut) {
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->Pindex() != trajParm_->Pindex()) return 0;

  // First frame setup - set up for the input parm, not necessarily the setup
  // parm; this allows things like atom strippping, etc. A stripped parm will
  // have the same pindex as the original parm.
  if (!trajIsOpen_) {
    if (debug_>0) rprintf("    Setting up %s for WRITE, %i atoms, originally %i atoms.\n",
                          trajName_,tparmIn->Natom(),trajParm_->Natom());
    trajParm_ = tparmIn;
    // Use parm to set up box info for the traj unless nobox was specified.
    // If box angles are present in traj they will be used instead.
    // NOTE: Probably not necessary to set box angles here, they are passed in
    if (!nobox_) {
      if (trajParm_->BoxType()!=Box::NOBOX) {
        trajio_->SetBox();
        //trajio->boxAngle[0]=trajParm->Box[3];
        //trajio->boxAngle[1]=trajParm->Box[4];
        //trajio->boxAngle[2]=trajParm->Box[5];
      }
    }
    // Set up write for the current parm file 
    if (trajio_->setupTrajout(trajParm_)) return 1;
    // Open output traj and mark as set up.
    if (trajio_->openTraj()) return 1;
    trajIsOpen_ = true;
    // If a framerange is defined set it to the begining of the range
    if (FrameRange_!=NULL) FrameRange_->Begin();
  }

  // If there is a framerange defined, check if this frame matches. If so,
  // write this frame and increment to the next frame in the range.
  if (FrameRange_!=NULL) {
    // If no more frames in the framerange, skip
    if ( FrameRange_->End() ) return 0;
    // Is this frame the next one in the range? 
    if ( FrameRange_->Current() != set ) return 0;
    FrameRange_->Next();
  }

  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (trajio_->writeFrame(set,FrameOut.X_,FrameOut.V_,FrameOut.box_,FrameOut.T_)) return 1;
  ++numFramesProcessed_;

  return 0;
}

// TrajectoryFile::TrajFilenameIs()
/** Call TrajectoryIO FilenameIs routine to check if input filename matches
  * full path of this trajectory file.
  */
bool TrajectoryFile::TrajFilenameIs(char *filenameIn) {
  return ( trajio_->FilenameIs(filenameIn) );
}

// TrajectoryFile::PrintInfo()
/** Print general trajectory information. Call trajio->Info for specific information.
  */
void TrajectoryFile::PrintInfo(int showExtended) {
  mprintf("  [%s] ",trajName_);
  trajio_->info();

  mprintf(", Parm %i",trajParm_->Pindex());

  if (trajio_->HasBox()) mprintf(" (with box info)");

  if (showExtended==0) {
    mprintf("\n");
    return;
  }

  if (fileAccess_==READTRAJ) {
    if (stop_!=-1 && total_frames_>0)
      //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
      mprintf(" (reading %i of %i)",total_read_frames_,total_frames_);
    else if (stop_!=-1 && total_frames_ < 0)
      mprintf(" (reading %i)",total_read_frames_);
    else
      mprintf(", unknown #frames, start=%i offset=%i",start_,offset_);
  } else {
    if (FrameRange_!=NULL)
      FrameRange_->PrintRange(": Writing frames",OUTPUTFRAMESHIFT);
    else
      mprintf(": Writing %i frames", trajParm_->Nframes());
    if (fileAccess_==APPENDTRAJ) mprintf(", appended"); 
  }
  if (debug_>0) 
    mprintf(", %i atoms, Box %i, seekable %i",trajParm_->Natom(),(int)trajio_->HasBox(),
            (int)trajio_->Seekable());
  mprintf("\n");
}

// Return private variables
int TrajectoryFile::CurrentFrame()       { return currentFrame_;       }
char *TrajectoryFile::TrajName()         { return (char*) trajName_;   }
Topology *TrajectoryFile::TrajParm()    { return trajParm_;           }
int TrajectoryFile::Start()              { return start_;              }
int TrajectoryFile::Total_Read_Frames()  { return total_read_frames_;  }
int TrajectoryFile::Total_Frames()       { return total_frames_;       }
int TrajectoryFile::NumFramesProcessed() { return numFramesProcessed_; }

// HasVelocity()
/** Return true if underlying trajio object indicates trajectory has velocity
  * information.
  */
bool TrajectoryFile::HasVelocity() { 
  if (trajio_!=NULL) return trajio_->HasVelocity(); 
  return false;
}

const char *TrajectoryFile::FileName() {
  if (trajio_!=NULL) return trajio_->Name();
  return 0;
}

