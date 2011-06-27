// RemdTraj
#include <cstdio> // sprintf
#include <cstring>
#include <cstdlib>
#include <cctype>
#include "RemdTraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
RemdTraj::RemdTraj() {
  Prefix=NULL;
  ExtWidth=0;
  CompressExt=NULL;
  repFilename=NULL;
  lowestRepnum=0;

  hasTrajout=false;
  remdtrajtemp=0.0;
  TemperatureList=NULL;
  // Used for writing replica trajectories
  remdX = NULL;
  remdbox[0]=0.0;
  remdbox[1]=0.0;
  remdbox[2]=0.0;
  remdbox[3]=0.0;
  remdbox[4]=0.0;
  remdbox[5]=0.0;
  remdT=0.0;
  remdN=0;
}

// DESTRUCTOR
RemdTraj::~RemdTraj() {
  std::vector<TrajectoryIO*>::iterator replica;
  //fprintf(stderr,"RemdTraj Destructor.\n");
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++)
    delete *replica;
  for (replica=REMDtrajout.begin(); replica!=REMDtrajout.end(); replica++)
    delete *replica;

  if (TemperatureList!=NULL) free(TemperatureList);
  if (Prefix!=NULL) free(Prefix);
  if (repFilename!=NULL) free(repFilename);
  if (CompressExt!=NULL) free(CompressExt);
}

/*
 * RemdTraj::SetupRead()
 * Set up this trajectory for processing as an REMD trajectory. All other
 * replica trajectory files are searched for, and the base filename, starting
 * replica, and number of replicas is recorded so the files may be opened
 * later during processing.
 */
/*
int RemdTraj::SetupRead(char *filenameIn, ArgList *argIn, AmberParm *parmIn) {
  int i,j, ExtWidth, startArg,stopArg,offsetArg;
  TrajectoryFile *traj;
  PtrajFile remdfile;  // Used to check for compression and extension
  int lastChar;        // The last character before compress extension (if any)
  char *ReplicaExt;    // Hold numeric extension
  char *CompressExt;   // Hold compression extension
  char *Prefix;        // Hold filename up until the numeric extension
  char *repFilename;   // Additional replica filenames
  char *lowestRepName; // Filename of the lowest # replica (passed in).

  // Check for lowest replica trajectory filename. Set to trajName
  if (filenameIn==NULL && argIn!=NULL)
    lowestRepName = argIn->getNextString();
  else
    lowestRepName = filenameIn;
  if (lowestRepName==NULL) {
    mprinterr("Error: RemdTraj::SetupRead: Filename is NULL.\n");
    return 1;
  }

  // Check for associated parm file
  if (parmIn==NULL) {
    mprinterr("Error: RemdTraj::SetupRead: Parm file is NULL.\n");
    return 1;
  }
  trajParm = parmIn;

  // Get target temperature
  remdtrajtemp=argIn->getKeyDouble("remdtrajtemp",0.0);

  // If remdout specified, treat all arguments following remdout as 
  // trajout arguments. Temperature trajectories will be written based 
  // on those arguments.
  if (argIn->hasKey("remdout")) {
    RemdOutArgs = argIn->SplitAt("remdout");
  }

  if (debug>0) {
    mprintf("    RemdTraj: Using specified file as lowest replica: %s\n",lowestRepName);
    mprintf("    RemdTraj: Frames at %lf K will be processed.\n",remdtrajtemp);
  }

  // Set up lowest replica file
  traj = new TrajectoryFile();
  if (traj==NULL) {
    mprinterr("Error: RemdTraj: Could not allocate memory for lowest replica.\n");
    return 1;
  }
  if (traj->SetupRead(lowestRepName,NULL,trajParm)) {
    mprinterr("    Error: RemdTraj: Could not set up lowest replica file %s\n",trajName);
    delete traj;
    return 1;
  }

  // Use this trajectory to set up overall stop, frames, and box info.
  // NOTE: Should check that this is the case for ALL frames.
  stop = traj->Total_Frames();
  total_frames = stop;
  //BoxType = T->BoxType;
  // Process start, stop, and offset args from user
  if (argIn!=NULL) {
    // Get any user-specified start, stop, and offset args
    // NOTE: For compatibility with ptraj start from 1
    startArg=argIn->getNextInteger(1);
    stopArg=argIn->getNextInteger(-1);
    offsetArg=argIn->getNextInteger(1);
    SetArgs(startArg,stopArg,offsetArg);
  }

  // Set replica traj name to base name of lowest replica
  SetTrajName( traj->TrajName() );
 
  // Check for temperature information
  if ( NoTempInfo(traj) ) return 1;

  // Add it to the list
  REMDtraj.push_back(traj);
  //REMDtraj.Info();
  numReplicas=1;

  // Scan for additional REMD traj files.
  // Assume the extension of given trajectory is the number of the lowest 
  // replica, and that the other files are in sequence (e.g. rem.000, rem.001, 
  // rem.002 or rem.000.gz, rem.001.gz, rem.002.gz etc).
  // Set up a PtrajFile here to check compression and extension.
  remdfile.SetupFile(lowestRepName, READ, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug);
  if (remdfile.compressType!=NONE && remdfile.Ext!=NULL) {
    CompressExt=remdfile.Ext; 
    lastChar = strlen(remdfile.filename) - strlen(remdfile.Ext);
  } else {
    CompressExt=(char*)"";
    lastChar = strlen(remdfile.filename);
  }

  // Find location of last '.' (not including compression extension) and store it in i
  i=-1;
  for (j=0; j<lastChar; j++)
    if (remdfile.filename[j]=='.') i=j;
  if (i==-1) {
    mprinterr("    Error: RemdTraj: Could not find numeric extension.\n");
    mprinterr("           Check that REMD files have naming scheme NAME.X\n");
    mprinterr("           where X is an integer of arbitrary width.\n");
    return 1;
  }

  // Store filename up until the numeric extension
  Prefix=(char*) malloc( (i+1) * sizeof(char));
  strncpy(Prefix,remdfile.filename,i);
  Prefix[i]='\0';
  //mprintf("  REMDDEBUG: Replica filename prefix: %s\n",Prefix);

  ExtWidth=lastChar - i - 1;
  //mprintf("  REMDDEBUG: Last . in %s located at %i\n",remdfile.filename,i);
  //mprintf("  REMDDEBUG: Allocating %i for extension\n",ExtWidth+1);
  //mprintf("  REMDDEBUG: EXTwidth=%i\n",ExtWidth);
  ReplicaExt=(char*) malloc( (ExtWidth+1) * sizeof(char));
  strncpy(ReplicaExt, remdfile.filename + i + 1, ExtWidth);
  ReplicaExt[ExtWidth]='\0'; 
  //mprintf("  REMDDEBUG: Replica extension is %s\n",ReplicaExt);

  // Check that all digits in extension are numbers 
  for (j=0; j<ExtWidth; j++) {
    if (isdigit(ReplicaExt[j])==0) {
      mprinterr("    Error: RemdTraj: Character #%i (%c) in extension %s is not a number!\n",
              j,ReplicaExt[j],ReplicaExt);
      free(ReplicaExt);
      free(Prefix);
      return 1;
    }
  }

  // Store lowest replica number
  j=atoi(ReplicaExt);
  //mprintf("  REMDDEBUG: index of first replica = %i\n",j);
  free(ReplicaExt);

  // Assume replica file names all have same length. 
  // Add some padding just in case
  repFilename=(char*) malloc( (strlen(remdfile.filename)+32) * sizeof(char));

  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j-1,CompressExt);
  if (fileExists(repFilename)) {
    mprintf("    Warning: RemdTraj: Replica# found lower than file specified with trajin!\n");
    mprintf("             (Found %s)\n",repFilename);
    mprintf("             trajin <file> remdtraj requires lowest # replica!\n");
  }   

  // Search for and add all replicas higher than this.
  j++;
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j,CompressExt);
  while ( fileExists(repFilename) ) {
    traj = new TrajectoryFile();
    if (traj==NULL) {
      mprinterr("Error: RemdTraj: Could not allocate memory for replica %0*i.\n",ExtWidth,j);
      return 1;
    }
    if (traj->SetupRead(lowestRepName,NULL,trajParm)) {
      mprinterr("    Error: RemdTraj: Could not set up replica file %s\n",repFilename);
      delete traj;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Check for temperature information
    if ( NoTempInfo(traj) ) {
      delete traj;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Check that #Frames and box info matches
    //if ( Frames!=T->Frames || BoxType!=T->BoxType ) {
    if ( total_frames!=traj->Total_Frames() ) {
      mprinterr("    Error: RemdTraj: #Frames (%i) in replica %i does not match\n",
                traj->Total_Frames(), j);
      mprinterr("           %i frames in lowest replica\n", total_frames);
      delete traj;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Set start/stop/offset args
    // Offset of 1 is for compatibility with ptraj - user frame #s are
    // offset +1 from internal cpptraj frame #s
    traj->SetArgs(start+1,stop,offset);
    // Add it to the list
    REMDtraj.push_back(traj);
    // Increment
    j++;
    numReplicas++;
    sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j,CompressExt);
  }
  free(repFilename);
  free(Prefix);
  return 0;
}
*/

/* RemdTraj::GetTemperatureList()
 * For each REMD input traj get the temperature of the first frame to
 * create a complete list of replica temperatures.
 */
int RemdTraj::SetupTemperatureList(int natom) {
  std::vector<TrajectoryIO*>::iterator replica;
  double *X, box[6];
  int repnum = 0;
  int err = 0;

  // Allocate temp space for reading in coords 
  X = (double*) malloc( natom * 3 * sizeof(double));
  if (X==NULL) return 1;
  // Allocate space for temperature list
  repnum = (int) REMDtraj.size();
  TemperatureList = (double*) malloc(repnum*sizeof(double));

  // Get a list of all temperatures present in input REMD trajectories
  repnum = 0;
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++) {
    err=0;
    if ( (*replica)->openTraj() ) {
      err = 1;
    } else {
      // Read 1 frame to get temperature
      err = ((*replica)->readFrame(0,X, box, TemperatureList+repnum));
      mprintf("      Rep %i T = %6.2lf\n",repnum,TemperatureList[repnum]);
      repnum++;
      // Close input traj
      (*replica)->closeTraj();
    }
    // If err!=0 an error occurred, bail out.
    if (err!=0) break;
  }
  free(X);
  return err;
}

/* RemdTraj::SetReplicaName()
 * Given the name of the lowest replica file, set basic replica name 
 * information. Get the file prefix up to the numerical extension, the 
 * width of the numerical extension, and the compression extension if 
 * the file is compressed.
 * Set and return the number of the lowest replica, or -1 on error.
 */
int RemdTraj::SetReplicaName(char* filenameIn) {
  PtrajFile remdfile;
  int lastChar,lastDot,j;
  char *ReplicaExt;

  // STEP 1 - Get filename Prefix
  // Assume the extension of given trajectory is the number of the lowest 
  // replica, and that the other files are in sequence (e.g. rem.000, rem.001, 
  // rem.002 or rem.000.gz, rem.001.gz, rem.002.gz etc).
  // Set up a PtrajFile here to check compression and extension.  
  if (remdfile.SetupFile(filenameIn, READ, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug)) 
    return -1;
  if (remdfile.compressType!=NONE && remdfile.Ext!=NULL) {
    CompressExt = (char*) malloc( (strlen(remdfile.Ext)+1) * sizeof(char));
    strcpy(CompressExt,remdfile.Ext);
    lastChar = strlen(remdfile.filename) - strlen(remdfile.Ext);
  } else {
    CompressExt = (char*) malloc( 1 * sizeof(char));
    strcpy(CompressExt,"");
    lastChar = strlen(remdfile.filename);
  }
  // Find location of last '.' (not including compression extension) and store it in i
  lastDot=-1;
  for (j=0; j<lastChar; j++)
    if (remdfile.filename[j]=='.') lastDot=j;
  if (lastDot==-1) {
    mprinterr("    Error: RemdTraj: Could not find numeric extension.\n");
    mprinterr("           Check that REMD files have naming scheme NAME.X\n");
    mprinterr("           where X is an integer of arbitrary width.\n");
    return -1;
  }
  // Store filename up until the numeric extension
  Prefix=(char*) malloc( (lastDot+1) * sizeof(char));
  strncpy(Prefix,remdfile.filename,lastDot);
  Prefix[lastDot]='\0';
  //mprintf("  REMDDEBUG: Replica filename prefix: %s\n",Prefix);

  // STEP 2 - Get the numerical extension
  ExtWidth=lastChar - lastDot - 1;
  //mprintf("  REMDDEBUG: Last . in %s located at %i\n",remdfile.filename,i);
  //mprintf("  REMDDEBUG: Allocating %i for extension\n",ExtWidth+1);
  //mprintf("  REMDDEBUG: EXTwidth=%i\n",ExtWidth);
  ReplicaExt=(char*) malloc( (ExtWidth+1) * sizeof(char));
  strncpy(ReplicaExt, remdfile.filename + lastDot + 1, ExtWidth);
  ReplicaExt[ExtWidth]='\0';
  //mprintf("  REMDDEBUG: Replica extension is %s\n",ReplicaExt);
  // Check that all digits in extension are numbers
  for (j=0; j<ExtWidth; j++) {
    if (isdigit(ReplicaExt[j])==0) {
      mprinterr("    Error: RemdTraj: Character #%i (%c) in extension %s is not a number!\n",
              j,ReplicaExt[j],ReplicaExt);
      free(ReplicaExt);
      return -1;
    }
  }
  // Store lowest replica number
  lowestRepnum=atoi(ReplicaExt);
  //mprintf("  REMDDEBUG: index of first replica = %i\n",j);
  free(ReplicaExt);

  // Assume replica file names all have same length. 
  // Add some padding just in case
  repFilename=(char*) malloc( (strlen(remdfile.filename)+32) * sizeof(char));

  return lowestRepnum;
}

/* RemdTraj::GetReplicaName()
 * If name information has already been set by SetRepicaName, return
 * the expected replica filename given a replica number.
 */
char *RemdTraj::GetReplicaName(int repnum) {
  if (repFilename==NULL) return NULL;
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,repnum,CompressExt);
  return repFilename;
}

/* RemdTraj::openTraj()
 * Open each trajectory in the list. 
 */
int RemdTraj::openTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // Open the trajectory
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++)
    (*replica)->openTraj();

  // Open output trajectories
  //if (RemdOutArgs!=NULL) {
  //  for (traj=REMDtrajout.begin(); traj!=REMDtrajout.end(); traj++)
  //    (*traj)->BeginTraj(false);
  //}

  return 0;
}

/* RemdTraj::closeTraj()
 * Close all trajectories in the list.
 */
void RemdTraj::closeTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // Close the trajectory
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++)
    (*replica)->closeTraj();

  // Close output trajectories
  //if (RemdOutArgs!=NULL) {
  //  for (traj=REMDtrajout.begin(); traj!=REMDtrajout.end(); traj++)
  //    (*traj)->EndTraj();
  //}
}

/* RemdTraj::readFrame()
 */
int RemdTraj::readFrame(int set, double *X, double *box, double *T) {
  std::vector<TrajectoryIO *>::iterator reptrajin;
  std::vector<TrajectoryIO *>::iterator reptrajout;
  int trajout, nrep;
  bool found = false;
  
  for (reptrajin=REMDtraj.begin(); reptrajin!=REMDtraj.end(); reptrajin++) {
    // No conversion to replica trajectories: Just find target temp
    if (!hasTrajout) {
      if ((*reptrajin)->readFrame(set,X,box,T)) return 1;
      // Check if this is the target temp
      if (*T == remdtrajtemp) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        return 0;
      }

    // All input REMD trajectories converted to temperature trajectories. 
    } else {
      // remdframe is allocated in SetupRead
      if ((*reptrajin)->readFrame(set,remdX,remdbox,&remdT)) return 1;
      // Check if this is the target temp. If so, set main Frame coords/box/temp
      if (remdT == remdtrajtemp) {
        //printf("REMDTRAJ: remdout: Set %i TEMP=%lf\n",set,remdframe->T);
        for (int x=0; x < remdN; x++)
          X[x] = remdX[x];
        for (int b=0; b < 6; b++)
          box[b] = remdbox[b];
        *T = remdT;
        found=true;
      }
      // Figure out which output file matches this temperature
      trajout=-1; // Will hold the index of target output traj for T
      nrep=0;
      for (reptrajout=REMDtrajout.begin(); reptrajout!=REMDtrajout.end(); reptrajout++) {
        if (TemperatureList[nrep] == remdT) {
          trajout = nrep;
          break;
        }
        nrep++;
      }
      if (trajout==-1) {
        mprinterr("\nError: RemdTraj: remdout: Temperature %6.2lf not found in Temperature list.\n",
                  remdT);
        return 1;
      }
      // Write to output traj
      (*reptrajout)->writeFrame(set, remdX,remdbox,remdT);
    }
  }  // END LOOP over input remd trajectories
  if (found) return 0;
  // If we have made it here this means target was not found
  mprinterr("\nError: RemdTraj: Final repTemp value read= %lf, set %i\n",*T,set);
  mprinterr("Could not find target %lf in any of the replica trajectories.\n",
            remdtrajtemp);
  mprinterr("Check that all replica trajectory files were found and that\n");
  mprinterr("none of the trajectories are corrupted (e.g. missing a temperature).\n");
  return 1;
}

/* RemdTraj::info()
 */
void RemdTraj::info() {
  mprintf("REMD trajectories (%i total, lowest replica: %s)\n",
          (int)REMDtraj.size(),GetReplicaName(lowestRepnum));
  if (hasTrajout) {
    mprintf("        remdout: trajectories will be converted to temperature trajectories:\n");
    //REMDtrajout.Info(8);
  }
  mprintf("        Looking for frames at %8.2lf K",remdtrajtemp);
}
    
