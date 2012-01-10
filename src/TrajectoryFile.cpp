// TrajectoryFile
#include <cstdlib> //div_t
#include <cstring> // TrajName
#include "TrajectoryFile.h"
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
#include "RemdTraj.h"
#include "Traj_CharmmDcd.h"

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
  trajectoryIsOpen=false;
  total_read_frames=-1;
  boxType=NOBOX;
  currentFrame=0;
  targetSet=0;
  frameskip=1;
  FrameRange=NULL;
  nobox=false;
}

// DESTRUCTOR
TrajectoryFile::~TrajectoryFile() {
  if (trajectoryIsOpen) this->EndTraj();
  if (trajio!=NULL) delete trajio;
  if (FrameRange!=NULL) delete FrameRange;
  if (progress!=NULL) delete progress;
  if (trajName!=NULL) free(trajName);
}

// TrajectoryFile::SetDebug()
/** Set debug level.  */
void TrajectoryFile::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("  TrajectoryFile debug level set to %i\n",debug);
}

// TrajectoryFile::SetTrajName()
/** Set trajName to be a copy of input string.  */
void TrajectoryFile::SetTrajName(char *nameIn) {
  trajName=(char*) realloc( trajName, (strlen(nameIn)+1) * sizeof(char) );
  strcpy(trajName,nameIn);
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
TrajectoryIO *TrajectoryFile::setupRemdTrajIO(char *lowestRepName, double remdtrajtemp, 
                                              char *remdout, FileFormat remdfmt) 
{
  RemdTraj *remdio=NULL;
  TrajectoryIO *replica0=NULL;
  char *repFilename;
  int repnum, repframes;

  // Add lowest replica to the list. Initial set up of remd trajio object
  remdio = new RemdTraj();
  remdio->remdtrajtemp = remdtrajtemp;
  remdio->seekable = trajio->seekable;
  remdio->hasBox = trajio->hasBox;
  remdio->hasVelocity = trajio->hasVelocity;
  remdio->hasTemperature = trajio->hasTemperature;
  remdio->boxAngle[0]=trajio->boxAngle[0];
  remdio->boxAngle[1]=trajio->boxAngle[1];
  remdio->boxAngle[2]=trajio->boxAngle[2];
  remdio->REMDtraj.push_back( trajio ); 

  // ------------------------------------------------------
  // Scan for additional REMD traj files.
  // Set base replica name information and get lowest replica number
  // MUST USE FULL PATH HERE, trajName is only the base filename
  repnum = remdio->SetReplicaName(lowestRepName);
  if (repnum < 0) {
    delete remdio;
    return NULL;
  }

  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  repFilename = remdio->GetReplicaName(repnum-1);
  if (fileExists(repFilename)) {
    mprintf("    Warning: RemdTraj: Replica# found lower than file specified with trajin!\n");
    mprintf("             (Found %s)\n",repFilename);
    mprintf("             trajin <file> remdtraj requires lowest # replica!\n");
  }

  // Search for and add all replicas higher than this.
  repnum++;
  repFilename = remdio->GetReplicaName(repnum);
  while ( fileExists(repFilename) ) {
    // Set up replica file repnum trajectory IO
    replica0 = setupTrajIO(repFilename,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE);
    if (replica0==NULL) {
      mprinterr("    Error: RemdTraj: Could not set up replica %i file %s\n",repnum,repFilename);
      delete remdio;
      return NULL;
    }
    // Pushing replica0 here allows the remdio destructor to handle it on errors
    remdio->REMDtraj.push_back( replica0 );
    // Check that number of frames matches
    repframes = replica0->setupRead(trajParm);
    if (repframes < 0 || repframes != total_frames) {
      mprinterr("    Error: RemdTraj: Replica %i frames (%i) does not match\n",
                repnum,repframes);
      mprinterr("           # frames in replica 0 (%i).\n",total_frames);
      delete remdio;
      return NULL;
    }
    // Check for temperature information
    if ( !replica0->hasTemperature) {
      mprinterr("    Error: RemdTraj: Replica %i does not have temperature info.\n",repnum);
      delete remdio;
      return NULL;
    }
    // Check box information
    if ( replica0->hasBox != trajio->hasBox ) {
      mprinterr("    Error: RemdTraj: Replica %i box info does not match lowest replica.\n");
      delete remdio;
      return NULL;
    }
    // Increment
    repnum++;
    repFilename = remdio->GetReplicaName(repnum);
  }

  // If remdout was specified, set up output trajectories
  if (remdout!=NULL) {
    // Set up temperature list
    if ( remdio->SetupTemperatureList(trajParm->natom) ) {
      mprinterr("Error: RemdTraj: remdout: could not get temperature list.\n");
      delete remdio;
      return NULL;
    }
    // BEGIN LOOP over Temperature
    for (repnum=0; repnum < (int)remdio->REMDtraj.size(); repnum++) { 
      // Set up output filename for this temperature
      repFilename = remdio->GetTemperatureName(remdout,repnum);
      mprintf("    Creating remd output traj: %s\n",repFilename);
      // Set up file with given type and format, and set up replica0 for the format.
      if ( (replica0=setupTrajIO(repFilename,WRITE,remdfmt,UNKNOWN_TYPE))==NULL ) {
        mprinterr("    Error: Could not set up T-replica file %s for writing.\n",repFilename);
        delete remdio;
        return NULL;
      }
      // Add to the remd trajout list
      remdio->REMDtrajout.push_back( replica0 );
      // Set up write here, trajParm will not change
      replica0->hasBox = trajio->hasBox;
      if (replica0->setupWrite(trajParm)) {
        delete remdio;
        return NULL;
      }
    } // END LOOP over input remd trajectories
  } // END REMDtrajout

  // Since the repeated calls to setupTrajIO have overwritten trajName, 
  // reset it to the lowest replica name. Note that unlike the non-replica
  // trajName, this will include the full path.
  SetTrajName( remdio->GetLowestReplicaName() );

  return (TrajectoryIO*) remdio;
}

// TrajectoryFile::setupTrajIO()
/** Set up basic trajectory file for given access type. Set the trajectory
  * name to be the base filename. Return the trajectory IO object for 
  * the format.
  */
TrajectoryIO *TrajectoryFile::setupTrajIO(char *tname, AccessType accIn,
                                        FileFormat fmtIn, FileType typeIn) {
  TrajectoryIO *tio = NULL;
  CpptrajFile *basicTraj = NULL;

  basicTraj = new CpptrajFile();
  if (basicTraj->SetupFile(tname,accIn,fmtIn,typeIn,debug)) {
    //mprinterr("    Error: Could not set up file %s.\n",tname);
    delete basicTraj;
    return NULL;
  }
  // Set trajectory name to the base filename
  SetTrajName( basicTraj->basefilename );

  // Allocate IO type based on format
  switch ( basicTraj->fileFormat ) {
    case AMBERRESTART: tio = new AmberRestart(); break;
    case AMBERTRAJ   : tio = new AmberCoord();    break;
    case AMBERNETCDF :
#ifdef BINTRAJ
      tio = new AmberNetcdf();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
#endif
      break;
    case AMBERRESTARTNC :
#ifdef BINTRAJ
      tio = new AmberRestartNC();
#else
      mprinterr("    Error: Can not set up trajectory (%s):\n",tname);
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
      delete basicTraj;
#endif
      break;
    case PDBFILE     : tio = new PDBfile();      break;
    case CONFLIB     : tio = new Conflib();      break;
    case MOL2FILE    : tio = new Mol2File();     break;
    case CHARMMDCD   : tio = new CharmmDcd();    break;
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",tname);
      delete basicTraj;
      return NULL;
  }

  // Happens when memory cannot be allocd, or not compiled for netcdf
  if (tio==NULL) return NULL;

  // Set debug level
  tio->SetDebug(debug);

  // Place the basic file in the trajectory IO class
  tio->SetFile(basicTraj);

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
  if (argIn==NULL) return 0;
  int startArg = argIn->getNextInteger(1);
  int stopArg = argIn->getNextInteger(-1);
  int offsetArg = argIn->getNextInteger(1);

#ifdef DEBUGTRAJ
  mprintf("DEBUG [%s] SetArgs: Original start, stop: %i %i\n",trajName,start,stop);
  mprintf("DEBUG [%s] SetArgs: Original startArg, stopArg: %i %i\n",trajName,startArg,stopArg);
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
    mprintf("DEBUG [%s] SetArgs: Start %i Stop %i  Offset %i\n",trajName,start,stop,offset);
  return 0;
}

// TrajectoryFile::SetBoxType()
/** Based on box angles in the given trajectory IO object and information in
  * the associated parmtop, set the box type for this trajectory. Return 0
  * if successful, 1 if an error occurs.
  */
// NOTE: Move to BoxType?
int TrajectoryFile::SetBoxType(TrajectoryIO *tio) {
  if (tio==NULL) return 1;
  if (trajParm==NULL) return 1;
  // If box coords present but no box info in associated parm, print
  // a warning.
  if (trajParm->boxType == NOBOX) {
    mprintf("\tWarning: Box info present in trajectory %s but not in\n",trajName);
    mprintf("\t         associated parm %s\n",trajParm->parmName);
  }
  boxType = CheckBoxType(tio->boxAngle,debug);
  // If box coords present but returned box type is NOBOX then no angles 
  // present in trajectory. Set box angles to parm default. If parm has
  // no box information then assume orthorhombic box.
  // NOTE: Is this only good for amber trajectory?
  if (boxType == NOBOX) {
    if (trajParm->boxType == NOBOX) {
      //mprinterr("Error: No angle information present in trajectory %s\n",trajName);
      //mprinterr("       or parm %s.\n",trajParm->parmName);
      //return 1;
      mprintf("\tWarning: No angle information present in trajectory %s\n",trajName);
      mprintf("\t         or parm %s - setting angles to 90.0!\n",trajParm->parmName);
      tio->boxAngle[0] = 90.0;
      tio->boxAngle[1] = 90.0;
      tio->boxAngle[2] = 90.0;
    } else {
      if (debug>0) {
        mprintf("\tWarning: No angle information present in trajectory %s:\n",trajName);
        mprintf("\t         Using angles from parm %s (beta=%lf).\n",trajParm->parmName,
                trajParm->Box[3]);
      }
      tio->boxAngle[0] = trajParm->Box[3];
      tio->boxAngle[1] = trajParm->Box[4];
      tio->boxAngle[2] = trajParm->Box[5];
    }
    // Set trajectory box type from angles in boxAngle
    boxType = CheckBoxType(tio->boxAngle,debug);
  }
  if (debug>0 || trajParm->boxType == NOBOX) {
    mprintf("\t[%s] Box type is",trajName);
    if (boxType==NOBOX) mprintf(" None.\n");
    else if (boxType==ORTHO) mprintf(" Orthorhombic.\n");
    else if (boxType==NONORTHO) mprintf(" NonOrthorhombic.\n");
  }
  // If no box info in parm, set it from trajectory
  if (trajParm->boxType == NOBOX) {
    mprintf("\tWarning: Setting parm %s box information from trajectory %s.\n",
            trajParm->parmName,trajName);
    trajParm->boxType = boxType;
    trajParm->Box[3] = tio->boxAngle[0]; 
    trajParm->Box[4] = tio->boxAngle[1]; 
    trajParm->Box[5] = tio->boxAngle[2]; 
  }
  return 0;
}

// TrajectoryFile::SingleFrame()
/** Tell the trajectory to set up stop and offset so that only start frame
  * will be processed.
  */
void TrajectoryFile::SingleFrame() {
  stop = start + 1;
  offset = 1;
  // Call setupFrameInfo to recalc total_read_frames. Since setupFrameInfo 
  // should have already been called in SetupRead (and thus any errors 
  // handled there) dont check for an error here. It should return 1.
  if ( setupFrameInfo() != 1 ) {
    mprintf("  Warning: Single frame requested for %s but not calcd!\n",trajName);
    mprintf("           start/stop/offset (%i, %i, %i)\n",start+1,stop+1,offset);
  }
}

// TrajectoryFile::SetupRead()
/** Set up trajectory for reading. Input trajectory filename can be specified
  * explicitly, or if not it should be the second argument in the given
  * argument list. Associate this trajectory with the given parm file.
  */
int TrajectoryFile::SetupRead(char *tnameIn, ArgList *argIn, AmberParm *tparmIn) {
  char *tname = NULL;
  // REMD
  double remdtrajtemp = 0.0;
  char *remdout = NULL;
  bool remdtraj=false;
  FileFormat remdfmt = AMBERTRAJ;

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
  trajParm = tparmIn;

  // Check for remdtraj keyword; if present, get args pertaining to REMD traj.
  if (argIn!=NULL && argIn->hasKey("remdtraj")) {
    remdtraj=true;
    // Get target temperature
    remdtrajtemp=argIn->getKeyDouble("remdtrajtemp",0.0);
    // If remdout specified, get a filename and optionally a format keyword
    // (default amber traj if none specified) that will be used for writing 
    // out temperature trajectories. Currently only AmberNetcdf and AmberTraj
    // formats supported for this option.
    remdout = argIn->getKeyString("remdout",NULL);
    if (remdout!=NULL) remdfmt=getFmtFromArg(argIn,AMBERTRAJ);
    if (remdfmt!=AMBERTRAJ && remdfmt!=AMBERNETCDF) {
      mprinterr("Error: remdout (%s): Unsupported format. Currently only amber\n",remdout);
      mprinterr("       trajectory and amber netcdf files supported for remdout.\n");
      return 1;
    }
  }

  // Set this trajectory access to READ
  fileAccess = READ;

  // Set up file; among other things this will determine the type and 
  // format, and set up trajio for the format.
  if ( (trajio=setupTrajIO(tname,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE))==NULL ) {
    mprinterr("    Error: Could not set up file %s for reading.\n",tname);
    return 1;
  }

  // Set up the format for reading and get the number of frames.
  // -1 indicates an error.
  // -2 indicates the number of frames could not be determined, read to EOF.
  total_frames = trajio->setupRead(trajParm);
  if (total_frames == -1) {
    mprinterr("    Error: Could not set up %s for reading.\n",tname);
    return 1;
  }
  if (total_frames>-1)
    mprintf("\t[%s] contains %i frames.\n",trajName,total_frames);
  else
    mprintf("\t[%s] contains an unknown number of frames.\n",trajName);

  // Set stop based on calcd number of frames.
  if (total_frames==0) {
    mprinterr("  Error: trajectory %s contains no frames.\n",trajName);
    return 1;
  }
  if (total_frames>0)
    stop = total_frames; 
  else
    stop = -1;

  // If the trajectory has box coords, set the box type from the box Angles.
  if (trajio->hasBox) {
    //mprintf("DEBUG:\tBOX ANGLES: %lf %lf %lf\n",trajio->boxAngle[0],
    //        trajio->boxAngle[1], trajio->boxAngle[2]);
    if (SetBoxType( trajio )) return 1;
  }

  // Set start, stop, and offset args from user input.
  if (SetArgs(argIn)) return 1;

  // Call setupFrameInfo to calc actual start and stop values based on
  // offset, as well as total_read_frames
  if ( setupFrameInfo() == 0 ) {
    mprinterr("  Error: No frames will be read from %s based on start, stop,\n",trajName);
    mprinterr("         and offset values (%i, %i, %i)\n",start+1,stop+1,offset);
    return 1;
  }

  // For replica trajectories, replace the current trajio (which should be
  // the lowest replica) with a special trajio object that will have 
  // a list of trajio objects, each one pertaining to a different replica.
  if (remdtraj) {
    // Check that this trajio has temperature info
    if (!trajio->hasTemperature) {
      mprinterr("Error: RemdTraj: Lowest replica file %s does not have temperature info.\n",
                trajName);
      return 1;
    }
    // Replace this trajio with a special replica traj one that will contain
    // trajio objects for all replicas
    if ( (trajio = setupRemdTrajIO(tname, remdtrajtemp, remdout, remdfmt))==NULL )
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
  div_t divresult;
  // DEBUG - No mpi for now
  int worldrank = 0;
  int worldsize = 1;

  //mprintf("DEBUG: Calling setupFrameInfo for %s with %i %i %i\n",trajName,
  //        start,stop,offset);
  if (stop==-1) return -1;

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

#ifdef TRAJDEBUG
  mprintf("DEBUG SETUPFRAMEINFO: %i-%i total %i\n",start,stop,total_read_frames);
#endif

  return total_read_frames;
}

// TrajectoryFile::getFmtFromArg()
/** Given an arglist, search for one of the file format keywords.
  * Default to def. 
  * NOTE: def should probably not be allowed to be UNKNOWN_FORMAT,
  * but this is currently not explicitly checked.
  */
FileFormat TrajectoryFile::getFmtFromArg(ArgList *argIn, FileFormat def) {
  FileFormat writeFormat = def;
  if (argIn==NULL) return writeFormat;
  if      ( argIn->hasKey("pdb")      ) writeFormat=PDBFILE;
  else if ( argIn->hasKey("data")     ) writeFormat=DATAFILE;
  else if ( argIn->hasKey("netcdf")   ) writeFormat=AMBERNETCDF;
  else if ( argIn->hasKey("restart")  ) writeFormat=AMBERRESTART;
  else if ( argIn->hasKey("ncrestart")) writeFormat=AMBERRESTARTNC;
  else if ( argIn->hasKey("restartnc")) writeFormat=AMBERRESTARTNC;
  else if ( argIn->hasKey("mol2")     ) writeFormat=MOL2FILE;
  else if ( argIn->hasKey("dcd")      ) writeFormat=CHARMMDCD;
  return writeFormat;
}

// TrajectoryFile::SetupWriteWithArgs()
/** Like SetupWrite, but intended for internal use. Allows a static
  * space-separated string to be passed in, which will be converted
  * to an argument list and passed to SetupWrite.
  */
int TrajectoryFile::SetupWriteWithArgs(char *tnameIn, const char *argstring,
                                       AmberParm *tparmIn, FileFormat fmtIn) {
  ArgList tempArg;
  char *tempString;
  //tempArg.SetDebug(2);
  // Since ArgList uses strtok cannot pass it a const string, copy to 
  // temporary string.
  tempString = new char[ strlen(argstring) + 1];
  strcpy(tempString,argstring);
  tempArg.SetList(tempString, " ");
  delete[] tempString;
  return SetupWrite(tnameIn,&tempArg,tparmIn,fmtIn);
}

// TrajectoryFile::SetupWrite()
/** Set up trajectory for writing. Output trajectory filename can be specified
  * explicitly, or if not it should be the second argument in the given
  * argument list. Associate with the given parm file initially, but setup for 
  * the format is performed on the first write call to accomodate state changes
  * like stripped atoms and so on. 
  * NOTE: make remdtraj a generic trigger for hasTemperature?
  */
int TrajectoryFile::SetupWrite(char *tnameIn, ArgList *argIn, AmberParm *tparmIn,
                               FileFormat writeFormatIn) {
  char *tname = NULL;
  AccessType access = WRITE;
  FileFormat writeFormat = writeFormatIn;
  FileType writeType = UNKNOWN_TYPE;
  char *onlyframes=NULL;

  // Mark as not yet set up
  trajectoryIsOpen=false;

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
  // Set the write file format from arglist; if no format keywords present
  // this will default to writeFormatIn
  writeFormat = getFmtFromArg(argIn, writeFormatIn);
  // If still unknown format default to amber trajectory
  if (writeFormat==UNKNOWN_FORMAT) writeFormat=AMBERTRAJ;

  // Check for append keyword
  if (argIn!=NULL && argIn->hasKey("append")) access = APPEND;

  // Set this trajectory access to the specified type
  fileAccess = access;

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

    // Process any write arguments specific to certain formats not related
    // to parm file. Options related to parm file are handled on the first
    // write in WriteFrame.
    if (trajio->processWriteArgs(argIn)) {
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
  if (fileAccess!=READ) return 0;

  // Open the trajectory
  if (trajio->openTraj()) {
    mprinterr("Error: TrajectoryFile::BeginTraj: Could not open %s\n",trajName);
    return 1;
  }
  trajectoryIsOpen=true;
  numFramesProcessed=0;

  // Set up a progress bar
  if (showProgress) progress = new ProgressBar(total_read_frames);

  // Determine what frames will be read
  targetSet=start;
  if (trajio->seekable) {
    frameskip = offset;
    currentFrame = start;
  } else {
    frameskip = 1;
    currentFrame = 0;
  }

  return 0;
}

// TrajectoryFile::PrintInfoLine()
/** Print a smaller amount of information than the PrintInfo function to
  * a single line. Give the trajectory name, start, stop, and offset.
  */
void TrajectoryFile::PrintInfoLine() {
  //rprintf( "----- [%s] (%i-%i, %i) -----\n",trajName,currentFrame+1,stop+1,offset);
  if (stop!=-1)
    rprintf( "----- [%s] (%i-%i, %i) -----\n",trajName,start+1,stop+1,offset);
  else
    rprintf( "----- [%s] (%i-EOF, %i) -----\n",trajName,start+1,offset);
}

// TrajectoryFile::EndTraj()
/** Close the trajectory. Should be called after all reads / writes completed.
  */
int TrajectoryFile::EndTraj() {
  trajio->closeTraj();
  trajectoryIsOpen=false;
  return 0;
}

// TrajectoryFile::GetNextFrame()
/** Get the next target frame from trajectory. Update the number of frames
  * read while getting to target (if traj is seekable this will always be 1).
  * \return 1 on successful read.
  * \return 0 if no frames could be read.
  */
int TrajectoryFile::GetNextFrame(Frame& FrameIn) { 
  bool tgtFrameFound;
#ifdef TRAJDEBUG
  mprinterr("Getting frame %i from %s (stop=%i)\n",currentFrame,trajName,stop);
#endif
  // If the current frame is out of range, exit
  if (currentFrame>stop && stop!=-1) return 0;
  //if (currentFrame>stop) return 0;
  if (progress!=NULL) 
    progress->Update(numFramesProcessed);
    //progress->PrintBar(currentFrame);

  tgtFrameFound=false;

  while ( !tgtFrameFound ) {
#ifdef TRAJDEBUG
    mprinterr("Attempting read of frame %i from %s\n",currentFrame,trajName);
#endif
    if (trajio->readFrame(currentFrame,FrameIn.X,FrameIn.V,FrameIn.box,&(FrameIn.T))) return 0;
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
int TrajectoryFile::WriteFrame(int set, AmberParm *tparmIn, Frame &FrameOut) {
  // Check that input parm matches setup parm - if not, skip
  if (tparmIn->pindex != trajParm->pindex) return 0;

  // First frame setup - set up for the input parm, not necessarily the setup
  // parm; this allows things like atom strippping, etc. A stripped parm will
  // have the same pindex as the original parm.
  if (!trajectoryIsOpen) {
    if (debug>0) rprintf("    Setting up %s for WRITE, %i atoms, originally %i atoms.\n",
                         trajName,tparmIn->natom,trajParm->natom);
    trajParm = tparmIn;
    // Use parm to set up box info for the traj unless nobox was specified.
    // If box angles are present in traj they will be used instead.
    // NOTE: Probably not necessary to set box angles here, they are passed in
    if (!nobox) {
      if (trajParm->boxType!=NOBOX) {
        trajio->hasBox=true;
        trajio->boxAngle[0]=trajParm->Box[3];
        trajio->boxAngle[1]=trajParm->Box[4];
        trajio->boxAngle[2]=trajParm->Box[5];
      }
    }
    // Set up write for the current parm file 
    if (trajio->setupWrite(trajParm)) return 1;
    // Open output traj and mark as set up.
    if (trajio->openTraj()) return 1;
    trajectoryIsOpen=true;
    // If a framerange is defined set it to the begining of the range
    if (FrameRange!=NULL) FrameRange->Begin();
  }

  // If there is a framerange defined, check if this frame matches. If so,
  // write this frame and increment to the next frame in the range.
  if (FrameRange!=NULL) {
    // If no more frames in the framerange, skip
    if ( FrameRange->End() ) return 0;
    // NOTE: For compatibility with ptraj user frame args start at 1
    if ( FrameRange->Current() - 1 != set ) return 0;
    FrameRange->Next();
  }

  // Write
  //fprintf(stdout,"DEBUG: %20s: Writing %i\n",trajName,set);
  if (trajio->writeFrame(set,FrameOut.X,FrameOut.V,FrameOut.box,FrameOut.T)) return 1;
  numFramesProcessed++;

  return 0;
}

// TrajectoryFile::TrajFilenameIs()
/** Call TrajectoryIO FilenameIs routine to check if input filename matches
  * full path of this trajectory file.
  */
bool TrajectoryFile::TrajFilenameIs(char *filenameIn) {
  return ( trajio->FilenameIs(filenameIn) );
}

// TrajectoryFile::PrintInfo()
/** Print general trajectory information. Call trajio->Info for specific information.
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
    if (stop!=-1 && total_frames>0)
      //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
      mprintf(" (reading %i of %i)",total_read_frames,total_frames);
    else if (stop!=-1 && total_frames < 0)
      mprintf(" (reading %i)",total_read_frames);
    else
      mprintf(", unknown #frames, start=%i offset=%i",start,offset);
  } else {
    mprintf(": Writing %i frames", trajParm->parmFrames);
    if (fileAccess==APPEND) mprintf(", appended"); 
  }
  if (debug>0) mprintf(", %i atoms, Box %i, seekable %i",trajParm->natom,boxType,trajio->seekable);
  mprintf("\n");
}

// HasVelocity()
/** Return true if underlying trajio object indicates trajectory has velocity
  * information.
  */
bool TrajectoryFile::HasVelocity() { 
  if (trajio!=NULL) return trajio->hasVelocity; 
  return false;
}
