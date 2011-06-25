// TrajectoryFile
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"
#include <cstdlib> //div_t
#include <cstring>
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
//#include "RemdTraj.h"
/*
#include "AmberTraj.h"
#ifdef BINTRAJ
#  include "AmberNetcdf.h"
#  include "AmberRestartNC.h"
#endif
#include "PDBfile.h"
#include "AmberRestart.h"
#include "Conflib.h"
#include "Mol2File.h"
*/

// CONSTRUCTOR
TrajectoryFile::TrajectoryFile() {
  debug = 0;
  progress=NULL;
  trajio=NULL;
  trajName=NULL;
  trajParm=NULL;
  fileAccess=READ;
  start=0;
  stop=-1;
  offset=1;
  total_frames=0;
  numFramesProcessed=0;
  total_read_frames=-1;
  boxType=NOBOX;
  currentFrame=0;
  targetSet=0;
  frameskip=1;
  FrameRange=NULL;
  nobox=false;
  setupForWrite=false;
};

// DESTRUCTOR
TrajectoryFile::~TrajectoryFile() {
  if (trajio!=NULL) delete trajio;
  if (FrameRange!=NULL) delete FrameRange;
  if (progress!=NULL) delete progress;
  if (trajName!=NULL) free(trajName);
}

/* TrajectoryFile::SetDebug()
 * Set debug level.
 */
void TrajectoryFile::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("  TrajectoryFile debug level set to %i\n",debug);
}

/* TrajectoryFile::SetTrajName()
 * Set trajName to be a copy of input string.
 */
void TrajectoryFile::SetTrajName(char *nameIn) {
  trajName=(char*) realloc( trajName, (strlen(nameIn)+1) * sizeof(char) );
  strcpy(trajName,nameIn);
}

/* TrajectoryFile::setupRemdTraj()
 * Set this trajectory up as an REMD trajectory. A special trajio object
 * will be set up that itself contains a list of trajio objects, each one
 * corresponding to a replica. Returns the trajio class of the lowest
 * replica.
 */
/*TrajectoryIO *TrajectoryFile::setupRemdTraj(char *lowestRepName, ArgList *argIn) {
  RemdTraj *remdio=NULL;
  TrajectoryIO *replica0=NULL;
  ArgList *RemdOutArgs=NULL;
  double remdtrajtemp;

  // Set up lowest replica file trajectory IO 
  replica0 = setupTrajIO(lowestRepName,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE);
  if (replica0==NULL) { 
    mprinterr("    Error: RemdTraj: Could not set up lowest replica file %s\n",lowestRepName);
    delete replica0;
    return NULL;
  }

  // Get target temperature
  remdtrajtemp=argIn->getKeyDouble("remdtrajtemp",0.0);

  // If remdout specified, treat all arguments following remdout as 
  // trajout arguments. Temperature trajectories will be written based 
  // on those arguments.
  if (argIn->hasKey("remdout")) {
    RemdOutArgs = argIn->SplitAt("remdout");
  }

  // Set up the lowest replica traj for reading. -1 indicates an error, 
  // -2 indicates # of frames could not be determined.
  total_frames = replica0->setupRead(trajParm->natom);
  if (total_frames == -1) {
    mprinterr("    Error: RemdTraj: Could not set up IO for %s\n",lowestRepName);
    delete replica0;
    return NULL;
  }

  // Ensure that lowest replica traj  has temperature information
  if (!replica0->hasTemperature) {
    mprinterr("    Error: RemdTraj: Lowest replica file %s does not have temperature info.\n",
              lowestRepName);
    delete replica0;
    return NULL;
  }

  // Set trajectory box type based on lowest replica traj box type
  if (SetBoxType(replica0)) {
    mprinterr("    Error: RemdTraj: Setting box info for lowest replica %s.\n",
              lowestRepName);
    delete replica0;
    return NULL;
  }

  // Process start, stop, and offset args
  if (SetArgs(argIn)) {
    delete replica0;
    return NULL;
  }

  // Add lowest replica to the list
  remdio = new RemdTraj();
  remdio->REMDtraj.push_back( replica0 ); 

  return (TrajectoryIO*) remdio;
}
*/

/* TrajectoryFile::setupTrajIO()
 * Set up basic trajectory file for read/write/append. Return the 
 * trajectory IO object for the format.
 */
TrajectoryIO *TrajectoryFile::setupTrajIO(char *tname, AccessType accIn,
                                        FileFormat fmtIn, FileType typeIn) {
  TrajectoryIO *tio = NULL;
  PtrajFile *basicTraj = NULL;

  fileAccess=accIn;
  basicTraj = new PtrajFile();
  if (basicTraj->SetupFile(tname,fileAccess,fmtIn,typeIn,debug)) {
    mprinterr("    Error: Could not set up file %s.\n",tname);
    delete basicTraj;
    return NULL;
  }

  // Allocate IO type based on format
  switch ( basicTraj->fileFormat ) {
//    case AMBERRESTART: tio = new AmberRestart(); break;
    case AMBERTRAJ   : tio = new AmberCoord();    break;
/*    case AMBERNETCDF :
#ifdef BINTRAJ
      tio = new AmberNetcdf();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      tio = new AmberRestartNC();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return NULL;
#endif
      break;
    case PDBFILE     : tio = new PDBfile();      break;
    case CONFLIB     : tio = new Conflib();      break;
    case MOL2FILE    : tio = new Mol2File();     break;
*/
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",tname);
      delete basicTraj;
      return NULL;
  }

  // Should only happen when memory cannot be allocd
  if (tio==NULL) return NULL;

  // Set debug level
  tio->SetDebug(debug);

  // Place the basic file in the trajectory IO class
  tio->SetFile(basicTraj);

  return tio;
}

/* TrajectoryFile::SetArgs()
 * Called after initial trajectory setup, set the start, stop, and offset args.
 * Do some bounds checking.
 * For compatibility with ptraj frames start at 1. So for a traj with 10 frames:
 * cpptraj: 0 1 2 3 4 5 6 7 8 9
 *   ptraj: 1 2 3 4 5 6 7 8 9 10
 * Defaults: startArg=1, stopArg=-1, offsetArg=1
 */
//void TrajectoryFile::SetArgs(int startArg, int stopArg, int offsetArg) {
int TrajectoryFile::SetArgs(ArgList *argIn) {
  // Set stop based on calcd number of Frames.
  if (total_frames == -2) {
    mprintf("  Warning: Could not predict # frames in %s. This usually indicates \n",trajName);
    mprintf("         a corrupted trajectory. Frames will be read until EOF.\n");
    stop=-1;
  } else if (total_frames==0) {
    mprinterr("  Error: trajectory %s contains no frames.\n",trajName);
    return 1;
  } else
    stop = total_frames;

  if (argIn==NULL) return 0;
  int startArg = argIn->getNextInteger(1);
  int stopArg = argIn->getNextInteger(-1);
  int offsetArg = argIn->getNextInteger(1);

#ifdef DEBUGTRAJ
  mprinterr("DEBUG: setArgs: Original start, stop: %i %i\n",startArg,stopArg);
#endif
  if (startArg!=1) {
    if (startArg<1) {
      mprintf("    Warning: %s start argument %i < 1, setting to 1.\n",trajName,startArg);
      start=0; // cpptraj = ptraj - 1
    } else if (total_frames>=0 && startArg>total_frames) {
      // If startArg==stopArg and is greater than # frames, assume we want
      // the last frame (useful when reading for reference structure).
      if (startArg==stopArg) {
        mprintf("    Warning: %s start %i > #Frames (%i), setting to last frame.\n",
                trajName,startArg,total_frames);
        start=total_frames - 1;
      } else {
        mprintf("    Warning: %s start %i > #Frames (%i), no frames will be processed.\n",
                trajName,startArg,total_frames);
        start=startArg - 1;
      }
    } else
      start=startArg - 1;
  }
  if (stopArg!=-1) {
    if ((stopArg - 1)<start) { // cpptraj = ptraj - 1
      mprintf("    Warning: %s stop %i < start, no frames will be processed.\n",
              trajName,stopArg);
      stop = start;
    } else if (total_frames>=0 && stopArg>total_frames) {
      mprintf("    Warning: %s stop %i >= #Frames (%i), setting to max.\n",
              trajName,stopArg,total_frames);
      stop=total_frames;
    } else
      stop=stopArg;
  }
  if (offsetArg!=1) {
    if (offset<1) {
      mprintf("    Warning: %s offset %i < 1, setting to 1.\n",
              trajName,offsetArg);
      offset=1;
    } else if (stop!=-1 && offsetArg > stop - start) {
      mprintf("    Warning: %s offset %i is so large that only 1 set will be processed.\n",
              trajName,offsetArg);
      offset=offsetArg;
    } else
      offset=offsetArg;
  }
  if (debug>0)
    mprintf("  [%s] Args: Start %i Stop %i  Offset %i\n",trajName,start,stop,offset);
  return 0;
}

/* TrajectoryFile::SetBoxType()
 * Based on box angles in the given trajectory IO object and information in
 * the associated parmtop, set the box type for this trajectory. Return 0
 * if successful, 1 if an error occurs.
 * NOTE: Move to BoxType?
 */
int TrajectoryFile::SetBoxType(TrajectoryIO *tio) {
  if (tio==NULL) return 1;
  if (trajParm==NULL) return 1;
  // If box coords present but no box info in associated parm, print
  // a warning.
  if (trajParm->boxType == NOBOX) {
    mprintf("Warning: Box info present in trajectory %s but not in\n",trajName);
    mprintf("         associated parm %s\n",trajParm->parmName);
  }
  boxType = CheckBoxType(tio->boxAngle,debug);
  // If box coords present but returned box type is NOBOX then no angles 
  // present in trajectory. Set box angles to parm default. If parm has
  // no box information then exit just to be safe.
  // NOTE: Is this only good for amber trajectory?
  if (boxType == NOBOX) {
    if (trajParm->boxType == NOBOX) {
      mprinterr("Error: No angle information present in trajectory %s\n",trajName);
      mprinterr("       or parm %s.\n",trajParm->parmName);
      return 1;
      //mprintf("         or parm %s - setting angles to 90.0!\n",trajParm->parmName);
      //tio->boxAngle[0] = 90.0;
      //tio->boxAngle[1] = 90.0;
      //tio->boxAngle[2] = 90.0;
    } else {
      if (debug>0) {
        mprintf("Warning: No angle information present in trajectory %s:\n",trajName);
        mprintf("         Using angles from parm %s (beta=%lf).\n",trajParm->parmName,
                trajParm->Box[3]);
      }
      tio->boxAngle[0] = trajParm->Box[3];
      tio->boxAngle[1] = trajParm->Box[4];
      tio->boxAngle[2] = trajParm->Box[5];
    }
    boxType = CheckBoxType(tio->boxAngle,debug);
  }
  if (debug>0) {
    mprintf("[%s] Box type is ",trajName);
    if (boxType==NOBOX) mprintf(" None.\n");
    else if (boxType==ORTHO) mprintf(" Orthorhombic.\n");
    else if (boxType==NONORTHO) mprintf(" NonOrthorhombic.\n");
  }
  return 0;
}

/* TrajectoryFile::SingleFrame()
 * Tell the trajectory to set up stop and offset so that only start frame
 * will be processed.
 */
void TrajectoryFile::SingleFrame() {
  stop = start;
  offset = 1;
}

/* TrajectoryFile::SetupRead()
 * Set up trajectory for reading. Input trajectory filename can be specified
 * explicitly, or if not it should be the second argument in the given
 * argument list. Associate this trajectory with the given parm file.
 */
int TrajectoryFile::SetupRead(char *tnameIn, ArgList *argIn, AmberParm *tparmIn) {
  char *tname = NULL;

  // Check for trajectory filename
  if (tnameIn==NULL && argIn!=NULL)
    tname = argIn->getNextString();
  else
    tname = tnameIn;
  if (tname==NULL) {
    mprinterr("Error: TrajectoryFile::SetupRead: Filename is NULL.\n");
    return 1;
  }
  // Set trajectory name 
  SetTrajName( tname );

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
  trajParm = tparmIn;

  // Check for remdtraj keyword; if present, set up as a Remd Trajectory.
  // This will set up a special trajio object. All further setup is performed
  // in setupRemdTraj.
//  if (argIn->hasKey("remdtraj")) {
//    if ( (trajio=setupRemdTraj(tname,argIn))==NULL ) {
//      mprinterr("    Error: Could not set up REMD trajectories using %s as lowest replica.\n",
//                tname);
//      return 1;
//    }
//    return 0;

  // Set up single file; among other things this will determine the type and 
  // format, and set up trajio for the format.
//  } else {
    if ( (trajio=setupTrajIO(tname,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE))==NULL ) {
      mprinterr("    Error: Could not set up file %s for reading.\n",tname);
      return 1;
    }
//  }

  // Set up the format for reading. -1 indicates an error, -2 indicates # of
  // frames could not be determined.
  total_frames = trajio->setupRead(trajParm->natom);
  if (total_frames == -1) {
    mprinterr("    Error: Could not set up IO for %s\n",tname);
    return 1;
  }

  // If the trajectory has box coords, set the box type from the box Angles.
  if (trajio->hasBox) {
    if (SetBoxType( trajio )) return 1;
  }

  // Set start, stop, and user args based on calcd number of frames
  if (SetArgs(argIn)) return 1;

  // Process any more arguments

  return 0;
}

/* TrajectoryFile::SetupFrameInfo()
 * Calculate number of frames that will be read based on start, stop, and
 * offset. Update the number of frames that will be read for the associated
 * traj parm. 
 * Return the total number of frames that will be read for this traj.
 */
int TrajectoryFile::SetupFrameInfo() {
  int Nframes;
  int ptraj_start_frame, ptraj_end_frame;
  int traj_start_frame, traj_end_frame;
  div_t divresult;
  // DEBUG - No mpi for now
  int worldrank = 0;
  int worldsize = 1;

  if (total_frames<=0) {
    //outputStart=-1;
    total_read_frames=0;
    return -1;
  }
  
  //mprintf("DEBUG: Calling setupFrameInfo for %s with %i %i %i\n",trajfilename,
  //        start,stop,offset);

  // Calc total frames that will be read
  // Round up
  divresult = div( (stop - start), offset);
  total_read_frames = divresult.quot;
  if (divresult.rem!=0) total_read_frames++;

  // Calc min num frames read by any given thread
  // last thread gets leftovers
  // In case of 0, last thread gets the frame
  divresult = div(total_read_frames,worldsize);
  Nframes=divresult.quot;

  // Ptraj (local) start and end frame
  ptraj_start_frame=(worldrank*Nframes);
  ptraj_end_frame=ptraj_start_frame+Nframes;
  // Last thread gets the leftovers
  if (worldrank==worldsize-1) ptraj_end_frame+=divresult.rem;

  // Actual Traj start and end frame (for seeking)
  traj_start_frame=(ptraj_start_frame*offset) + start;
  traj_end_frame=((ptraj_end_frame-1)*offset) + start;

  start=traj_start_frame;
  stop=traj_end_frame;

  trajParm->parmFrames+=total_read_frames;

  return total_read_frames;
}
  

/* TrajectoryFile::SetupWrite()
 * Set up trajectory for writing. Output trajectory filename can be specified
 * explicitly, or if not it should be the second argument in the given
 * argument list. Associate with the given parm file initially, but setup for 
 * the format is performed on the first write call to accomodate state changes
 * like stripped atoms and so on. 
 */
int TrajectoryFile::SetupWrite(char *tnameIn, ArgList *argIn, AmberParm *tparmIn) {
  char *tname = NULL;
  AccessType access = WRITE;
  FileFormat writeFormat = AMBERTRAJ;
  FileType writeType = UNKNOWN_TYPE;
  char *onlyframes=NULL;

  // Check for trajectory filename
  if (tnameIn==NULL && argIn!=NULL)
    tname = argIn->getNextString();
  else
    tname = tnameIn;
  if (tname==NULL) {
    mprinterr("Error: TrajectoryFile::SetupWrite: Filename is NULL.\n");
    return 1;
  }
  // Set trajectory name 
  SetTrajName( tname );

  // Check for associated parm file
  if (tparmIn==NULL) {
    mprinterr("Error: TrajectoryFile::SetupWrite: Parm file is NULL.\n");
    return 1;
  }
  trajParm = tparmIn;

  // Process arguments related to access and format
  if (argIn!=NULL) {
    // Check for append keyword
    if ( argIn->hasKey("append") ) access = APPEND;

    // Set the write file format
    if      ( argIn->hasKey("pdb")      ) writeFormat=PDBFILE;
    else if ( argIn->hasKey("data")     ) writeFormat=DATAFILE;
    else if ( argIn->hasKey("netcdf")   ) writeFormat=AMBERNETCDF;
    else if ( argIn->hasKey("restart")  ) writeFormat=AMBERRESTART;
    else if ( argIn->hasKey("ncrestart")) writeFormat=AMBERRESTARTNC;
    else if ( argIn->hasKey("mol2")     ) writeFormat=MOL2FILE;
  }

  // Set up file with given type and format, and set up trajio for the format.
  if( (trajio=setupTrajIO(tname,access,writeFormat,writeType))==NULL ) {
    mprinterr("    Error: Could not set up file %s for writing.\n",tname);
    return 1;
  }

  // Process additional arguments
  if (argIn!=NULL) {
    // Get specified title if any - will not set if NULL
    trajio->SetTitle( argIn->getKeyString("title", NULL) );

    // Get a frame range for trajout
    onlyframes = argIn->getKeyString("onlyframes",NULL);
    if (onlyframes!=NULL) {
      FrameRange = new Range();
      if ( FrameRange->SetRange(onlyframes) ) {
        mprintf("Warning: trajout %s: onlyframes: %s is not a valid range.\n",tname,onlyframes);
        delete FrameRange;
      } else {
        FrameRange->PrintRange("      Saving frames",0);
      }
    }

    // Check for nobox argument - will override any box info present in parm
    // when trajectory IO is set up.
    nobox = argIn->hasKey("nobox");

    // Process any write arguments specific to certain formats
    if (trajio->TrajFormat() == AMBERTRAJ) { 
      AmberCoord *T0 = (AmberCoord*) trajio;
      if (argIn->hasKey("remdtraj")) T0->SetRemdTraj();
    }
  }

  // No more setup here; Write is set up when first frame written.
  return 0;
}

/* TrajectoryFile::BeginTraj()
 * Prepare trajectory file for reading/writing. For reads, open the traj and 
 * set up a progress bar if requested. For writes, just open the file.
 */
int TrajectoryFile::BeginTraj(bool showProgress) {
  // Open the trajectory
  if (trajio->openTraj()) {
    mprinterr("Error: TrajectoryFile::BeginTraj: Could not open %s\n",trajName);
    return 1;
  }
  numFramesProcessed=0;

  // If writing, this is all that is needed
  if (fileAccess!=READ) return 0;

  // Set up a progress bar
  if (showProgress) progress = new ProgressBar(stop);

  // Determine what frames will be read
  targetSet=start;
  if (trajio->seekable) {
    frameskip = offset;
    currentFrame = start;
  } else {
    frameskip = 1;
    currentFrame = 0;
  }
  rprintf( "----- [%s] (%i-%i, %i) -----\n",trajName,currentFrame+1,stop+1,offset);

  return 0;
}

/* TrajectoryFile::EndTraj()
 * Close the trajectory. 
 */
int TrajectoryFile::EndTraj() {
  trajio->closeTraj();
  return 0;
}

/* TrajectoryFile::GetNextFrame()
 * Get the next target frame from trajectory. Update the number of frames
 * read while getting to target (if traj is seekable this will always be 1).
 * Return 1 on successful read.
 * Return 0 if no frames could be read.
 */
int TrajectoryFile::GetNextFrame(double *X, double *box, double *T) { 
  bool tgtFrameFound;
#ifdef TRAJDEBUG
  mprinterr("Getting frame %i from %s (stop=%i)\n",currentFrame,trajName,stop);
#endif
  // If the current frame is out of range, exit
  if (currentFrame>stop && stop!=-1) return 0;
  if (progress!=NULL) 
    progress->Update(currentFrame);
    //progress->PrintBar(currentFrame);

  tgtFrameFound=false;

  while ( !tgtFrameFound ) {
#ifdef TRAJDEBUG
    mprinterr("Attempting read of frame %i from %s\n",currentFrame,trajName);
#endif
    if (trajio->readFrame(currentFrame,X,box,T)) return 0;
    //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
#ifdef TRAJDEBUG
    mprinterr("Frame %i has been read from %s (target=%i)\n",currentFrame,trajName,targetSet);
#endif
    if (currentFrame==targetSet) {
      tgtFrameFound=true;
      targetSet+=offset;
    }
    numFramesProcessed++;
    currentFrame+=frameskip;
  }

  return 1;
}

/* TrajectoryFile::WriteFrame()
 * Write the given coordinates, box, and Temperature. If no parm present
 * set up trajectory to write for this parm. If trajectory has been
 * set up, only write if the given parm matches the setup parm.
 */
int TrajectoryFile::WriteFrame(int set, AmberParm *tparmIn, double *X,
                               double *box, double T) {
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->pindex != trajParm->pindex) return 0;

  // First frame setup - set up for the input parm, not necessarily the setup
  // parm; this allows things like atom strippping, etc. A stripped parm will
  // have the same pindex as the original parm.
  if (!setupForWrite) {
    if (debug>0) rprintf("    Setting up %s for WRITE, %i atoms, originally %i atoms.\n",
                         trajName,tparmIn->natom,trajParm->natom);
    trajParm = tparmIn;
    // Use parm to set up box info unless nobox was specified.
    if (!nobox) {
      if (trajParm->boxType!=NOBOX) {
        trajio->hasBox=true;
        trajio->boxAngle[0]=trajParm->Box[4];
        trajio->boxAngle[1]=trajParm->Box[5];
        trajio->boxAngle[2]=trajParm->Box[6];
      }
    } 
    if (trajio->setupWrite(trajParm->natom)) return 1;
    if (trajio->openTraj()) return 1;
    setupForWrite=true;
  }

  // If there is a framerange defined, check if this frame matches. If so, pop
  if (FrameRange!=NULL) {
    // If no more frames in the framerange, skip
    if ( FrameRange->empty() ) return 0;
    // NOTE: For compatibility with ptraj user frame args start at 1
    if ( FrameRange->front() - 1 != set ) return 0;
    FrameRange->pop_front();
  }

  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (trajio->writeFrame(set,X,box,T)) return 1;
  numFramesProcessed++;

  return 0;
}

/* TrajectoryFile::TrajFilenameIs()
 * Call TrajectoryIO FilenameIs routine to check if input filename matches
 * full path of this trajectory file.
 */
bool TrajectoryFile::TrajFilenameIs(char *filenameIn) {
  return ( trajio->FilenameIs(filenameIn) );
}

/* TrajectoryFile::PrintInfo()
 * Print general trajectory information. Call trajio->Info for specific information.
 */
void TrajectoryFile::PrintInfo(int showExtended) {
  mprintf("  [%s] ",trajName);
  trajio->info();

  mprintf(", Parm %i",trajParm->pindex);

  if (trajio->hasBox) mprintf(" (with box info)");

  if (showExtended==0) {
    mprintf("\n");
    return;
  }

  if (fileAccess==READ) {
    if (stop!=-1)
      //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
      mprintf(" (reading %i of %i)",total_read_frames,total_frames);
    else
      mprintf(", unknown #frames, start=%i offset=%i",start,offset);
  } else {
    mprintf(": Writing %i frames", trajParm->parmFrames);
    if (fileAccess==APPEND) mprintf(", appended"); 
  }
  if (debug>0) mprintf(", %i atoms, Box %i, seekable %i",trajParm->natom,boxType,trajio->seekable);
  mprintf("\n");
}

