// TrajectoryFile
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"
#include <cstdlib> //div_t
#include <cstring>
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
/*
#include "AmberTraj.h"
#ifdef BINTRAJ
#  include "AmberNetcdf.h"
#  include "AmberRestartNC.h"
#endif
#include "PDBfile.h"
#include "AmberRestart.h"
#include "RemdTraj.h"
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
  start=1;
  stop=-1;
  offset=1;
  total_frames=0;
  numFramesRead=0;
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

/* TrajectoryFile::HasTemperature
 * True if trajio indicates trajectory has temperature information.
 */
bool TrajectoryFile::HasTemperature() {
  if (trajio==NULL) return false;
  return trajio->hasTemperature;
}

/* TrajectoryFile::setupTraj()
 * Set up basic trajectory file for read/write/append. Set up the 
 * trajectory IO object for the format.
 */
int TrajectoryFile::setupTraj(char *tname, AccessType accIn,
                                  FileFormat fmtIn, FileType typeIn) {
  PtrajFile *basicTraj = NULL;

  fileAccess=accIn;
  basicTraj = new PtrajFile();
  if (basicTraj->SetupFile(tname,fileAccess,fmtIn,typeIn,debug)) {
    delete basicTraj;
    if (debug>1)
      mprinterr("    Error: Could not set up file %s.\n",tname);
    return 1;
  }

  // Set trajectory name to the base filename
  SetTrajName( basicTraj->basefilename );

  // Allocate IO type based on format
  switch ( basicTraj->fileFormat ) {
//    case AMBERRESTART: trajio = new AmberRestart(); break;
    case AMBERTRAJ   : trajio = new AmberCoord();    break;
/*    case AMBERNETCDF :
#ifdef BINTRAJ
      trajio = new AmberNetcdf();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return 1;
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      trajio = new AmberRestartNC();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
      return 1;
#endif
      break;
    case PDBFILE     : trajio = new PDBfile();      break;
    case CONFLIB     : trajio = new Conflib();      break;
    case MOL2FILE    : trajio = new Mol2File();     break;
*/
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",tname);
      delete basicTraj;
      return 1;
  }

  // Should only happen when memory cannot be allocd
  if (trajio==NULL) return 1;

  // Place the basic file in the trajectory IO class
  trajio->SetFile(basicTraj);

  return 0;
}

/* TrajectoryFile::SetArgs()
 * Called after initial trajectory setup, set the start, stop, and offset args.
 * Do some bounds checking.
 * For compatibility with ptraj frames start at 1. So for a traj with 10 frames:
 * cpptraj: 0 1 2 3 4 5 6 7 8 9
 *   ptraj: 1 2 3 4 5 6 7 8 9 10
 * Defaults: startArg=1, stopArg=-1, offsetArg=1
 */
void TrajectoryFile::SetArgs(int startArg, int stopArg, int offsetArg) {
  //mprintf("DEBUG: setArgs: Original start, stop: %i %i\n",startArg,stopArg);
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
  int startArg, stopArg, offsetArg;

  // Check for trajectory filename
  if (tnameIn==NULL && argIn!=NULL)
    tname = argIn->getNextString();
  else
    tname = tnameIn;
  if (tname==NULL) {
    mprinterr("Error: TrajectoryFile::SetupRead: Filename is NULL.\n");
    return 1;
  }

  // Check for associated parm file
  if (tparmIn==NULL) {
    mprinterr("Error: TrajectoryFile::SetupRead: Parm file is NULL.\n");
    return 1;
  }
  trajParm = tparmIn;

  // Check for remdtraj keyword; if present, set up as a Remd Trajectory.
  // This will set up a special trajio object.
  if (argIn->hasKey("remdtraj")) {
    if (setupRemdTraj(tname,argIn)) {
      mprinterr("    Error: Could not set up REMD trajectories using %s as lowest replica.\n",
                tname);
      return 1;
    }

  // Set up single file; among other things this will determine the type and 
  // format, and set up trajio for the format.
  } else {
    if (setupTraj(tname,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE)) {
      mprinterr("    Error: Could not set up file %s for reading.\n",tname);
      return 1;
    }
  }

  // Set up the format for reading. -1 indicates an error, -2 indicates # of
  // frames could not be determined.
  total_frames = trajio->setupRead(trajParm->natom);
  if (total_frames == -1) {
    mprinterr("    Error: Could not set up IO for %s\n",tname);
    return 1;
  }

  // If the trajectory has box coords, set the box type from the box Angles.
  if (trajio->hasBox) {
    // If box coords present but no box info in associated parm, print
    // a warning.
    if (trajParm->boxType == NOBOX) {
      mprintf("Warning: Box info present in trajectory %s but not in\n",tname);
      mprintf("         associated parm %s\n",trajParm->parmName);
    }
    boxType = CheckBoxType(trajio->boxAngle,debug);
    // If box coords present but returned box type is NOBOX then no angles 
    // present in trajectory. Set box angles to parm default. If parm has
    // no box information then exit just to be safe.
    // NOTE: Is this only good for amber trajectory?
    if (boxType == NOBOX) {
      if (trajParm->boxType == NOBOX) {
        mprinterr("Error: No angle information present in trajectory %s\n",tname);
        mprinterr("       or parm %s.\n",trajParm->parmName);
        return 1;
        //mprintf("         or parm %s - setting angles to 90.0!\n",trajParm->parmName);
        //trajio->boxAngle[0] = 90.0;
        //trajio->boxAngle[1] = 90.0;
        //trajio->boxAngle[2] = 90.0;
      } else {
        if (debug>0) {
          mprintf("Warning: No angle information present in trajectory %s:\n",tname);
          mprintf("         Using angles from parm %s (beta=%lf).\n",trajParm->parmName,
                  trajParm->Box[3]);
        }
        trajio->boxAngle[0] = trajParm->Box[3];
        trajio->boxAngle[1] = trajParm->Box[4];
        trajio->boxAngle[2] = trajParm->Box[5];
      }
      boxType = CheckBoxType(trajio->boxAngle,debug);
    }
  }

  // Set stop based on calcd number of Frames.
  if (total_frames == -2) {
    mprintf("  Warning: Could not predict # frames in %s. This usually indicates \n",tname);
    mprintf("         a corrupted trajectory. Frames will be read until EOF.\n");
    stop=-1;
  } else if (total_frames==0) {
    mprinterr("  Error: trajectory %s contains no frames.\n",tname);
    return 1;
  } else
    stop = total_frames;

  // Process any more arguments
  if (argIn!=NULL) {
    // Get any user-specified start, stop, and offset args
    // NOTE: For compatibility with ptraj start from 1
    startArg=argIn->getNextInteger(1);
    stopArg=argIn->getNextInteger(-1);
    offsetArg=argIn->getNextInteger(1);
    SetArgs(startArg,stopArg,offsetArg);

    // Process arguments related to specific formats
  }

  return 0;
}

/* TrajectoryFile::SetupFrameInfo()
 * Calculate number of frames that will be read based on start, stop, and
 * offset. Update the number of frames that will be read for the associated
 * traj parm. 
 * Return the total number of frames that will be read for this traj.
 */
int TrajectoryFile::SetupFrameInfo() {
  div_t divresult;

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
  if( setupTraj(tname,access,writeFormat,writeType) ) {
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
  numFramesRead=0;

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
  // If the current frame is out of range, exit
  if (currentFrame>stop && stop!=-1) return 0;
  if (progress!=NULL) 
    progress->Update(currentFrame);
    //progress->PrintBar(currentFrame);

  tgtFrameFound=false;

  while ( !tgtFrameFound ) {
    if (trajio->readFrame(currentFrame,X,box,T)) return 0;
    //printf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
#ifdef DEBUG
    dbgprintf("DEBUG:\t%s:  current=%i  target=%i\n",trajName,currentFrame,targetSet);
#endif
    if (currentFrame==targetSet) {
      tgtFrameFound=true;
      targetSet+=offset;
    }
    numFramesRead++;
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

