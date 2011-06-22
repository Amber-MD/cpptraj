// TrajectoryFile
#include "TrajectoryFile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
TrajectoryFile::TrajectoryFile() {
  debug = 0;
  trajio=NULL;
  progress=NULL;
  trajName=NULL;
  trajParm=NULL;
  fileAccess=READ;
  start=1;
  stop=-1;
  offset=1;
  Frames=0;
  boxType=NONE;
  FrameRange=NULL;
  nobox=false;
};

// DESTRUCTOR
TrajectoryFile::~TrajectoryFile() {
  if (trajio!=NULL) delete trajio;
  if (progress!=NULL) delete progress;
  //if (trajName!=NULL) free(trajName);
  if (FrameRange!=NULL) delete FrameRange;
}

/* TrajectoryFile::setupTraj()
 * Set up basic trajectory file for read/write/append. Set the trajectory IO 
 * object for the type.
 */
int TrajectoryFile::setupTraj(char *tname, AccessType accIn, FileFormat fmtIn, FileType typeIn) {
  PtrajFile *basicTraj;

  fileAccess=accIn;
  basicTraj = new PtrajFile(); 
  if (basicTraj->SetupFile(tname,fileAccess,fmtIn,typeIn,debug)) {
    delete basicTraj;
    if (debug>1)
      mprinterr("    Error: Could not set up file %s.\n",tname);
    return 1;
  }

  // Set trajectory name to the base filename
  trajName = basicTraj->basefilename;

/*
  // Allocate IO type based on format
  switch ( basicTraj->fileFormat ) {
    case AMBERRESTART: trajio = new AmberRestart(); break;
    case AMBERTRAJ   : trajio = new AmberTraj();    break;
    case AMBERNETCDF :
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
    default:
      mprinterr("    Error: Could not determine trajectory file %s type, skipping.\n",tname);
      delete basicTraj;
      return 1;
  }
*/
  
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
  //mprintf("DEBUG: SetArgs: Original start, stop: %i %i\n",startArg,stopArg);
  if (startArg!=1) {
    if (startArg<1) {
      mprintf("    Warning: %s start argument %i < 1, setting to 1.\n",trajName,startArg);
      start=0; // cpptraj = ptraj - 1
    } else if (Frames>=0 && startArg>Frames) {
      // If startArg==stopArg and is greater than # frames, assume we want
      // the last frame (useful when reading for reference structure).
      if (startArg==stopArg) {
        mprintf("    Warning: %s start %i > #Frames (%i), setting to last frame.\n",
                trajName,startArg,Frames);
        start=Frames - 1;
      } else {
        mprintf("    Warning: %s start %i > #Frames (%i), no frames will be processed.\n",
                trajName,startArg,Frames);
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
    } else if (Frames>=0 && stopArg>Frames) {
      mprintf("    Warning: %s stop %i >= #Frames (%i), setting to max.\n",
              trajName,stopArg,Frames);
      stop=Frames;
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

  // Set up file; among other things this will determine the type and format,
  // and set up trajio for the format.
  if (setupTraj(tname,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE)) {
    mprinterr("    Error: Could not set up file %s for reading.\n",tname);
    return 1;
  }

  // Set up the format for reading.
  if (trajio->setupRead(trajParm->natom)) {
    mprinterr("    Error: Could not set up IO for %s\n",tname);
    return 1;
  }

  if (argIn!=NULL) {
    // Get any user-specified start, stop, and offset args
    // NOTE: For compatibility with ptraj start from 1
    startArg=argIn->getNextInteger(1);
    stopArg=argIn->getNextInteger(-1);
    offsetArg=argIn->getNextInteger(1);
    SetArgs(startArg,stopArg,offsetArg);
  }

  return 0;
}

/* TrajectoryFile::SetupWrite()
 * Set up trajectory for writing. Output trajectory filename can be specified
 * explicitly, or if not it should be the second argument in the given
 * argument list. Setup for the format is not performed here, but on the 
 * first write call. 
 */
int TrajectoryFile::SetupWrite(char *tnameIn, ArgList *argIn) {
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
  if (setupTraj(tname,access,writeFormat,writeType)) {
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
    //T->WriteArgs(A);
  }

  // No more setup here; Write is set up when first frame written.
  return 0;
}

int BeginTraj() {return 1;}
int EndTraj() {return 1;}
int GetNextFrame() { return 1;}
int WriteFrame(int set, AmberParm *tparmIn) {return 1;}
