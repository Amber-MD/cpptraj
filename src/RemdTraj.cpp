// RemdTraj
#include <cstdio> // sprintf
#include <cstring>
#include <cstdlib>
#include <cctype>
#include "RemdTraj.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
RemdTraj::RemdTraj() {
  remdtrajtemp=0.0;
  numReplicas=0;
  RemdOutArgs=NULL;
  TemperatureList=NULL;
  remdframe=NULL;
}

// DESTRUCTOR
RemdTraj::~RemdTraj() {
  //fprintf(stderr,"RemdTraj Destructor.\n");
  if (RemdOutArgs!=NULL) delete RemdOutArgs;
  if (TemperatureList!=NULL) free(TemperatureList);
  if (remdframe!=NULL) delete remdframe;
}

/* RemdTraj::NoTempInfo()
 * Check that trajectory has temperature info.
 * Return 1 if no temperature info, 0 otherwise.
 */
bool RemdTraj::NoTempInfo(TrajectoryFile *T) {
  if (T==NULL) return true;
  if (!T->HasTemperature()) {
    mprinterr("Error: REMDTRAJ: Trajectory %s does not contain temperature information.\n",
            T->TrajName());
    return true;
  }
  return false;
}

/*
 * RemdTraj::SetupRead()
 * Set up this trajectory for processing as an REMD trajectory. All other
 * replica trajectory files are searched for, and the base filename, starting
 * replica, and number of replicas is recorded so the files may be opened
 * later during processing.
 */
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

  // If remdout was specified, set up output trajectories
  if (RemdOutArgs!=NULL) {
    // Use Prefix to hold copy of base filename, which will be replaced by filename.temperature
    Prefix = RemdOutArgs->CopyArg(1);
    // Allocate space for output filename. Should be enough room for
    // filename.temperature, allocate some extra space for future-proofing.
    repFilename=(char*) malloc( (strlen(Prefix)+32) * sizeof(char));
    // Allocate frame for output trajectories
    remdframe = new Frame(trajParm->natom, NULL);
    // Allocate space for temperature list
    TemperatureList = (double*) malloc(numReplicas*sizeof(double));
    j=0;
    // BEGIN LOOP over input remd trajectories
    // Get a list of all temperatures present in input REMDtraj
    for (std::list<TrajectoryFile *>::iterator it=REMDtraj.begin(); it!=REMDtraj.end(); it++) {
      (*it)->BeginTraj(false);
      // Read 1 frame to get temperature
      (*it)->GetNextFrame(remdframe->X, remdframe->box, &(remdframe->T));
      TemperatureList[j++] = remdframe->T;
      if (debug>0) mprintf("    Rep %i T = %6.2lf\n",j-1,TemperatureList[j-1]);
      // Set up output filename for this temperature
      sprintf(repFilename,"%s.%6.2lf",Prefix,remdframe->T);
      //mprintf("    Creating remd output traj: %s\n",repFilename);
      // Replace the filename arg at position 1 with filename.Temperature
      if ( RemdOutArgs->ReplaceArg(1,repFilename) ) {
        mprinterr("Error: RemdTraj::setup(): Could not set output replica traj filename [%s].\n",
                  repFilename);
        j=-1;
      }
      //RemdOutArgs->print(); // DEBUG
      // Set up output remd traj file. Reset the arglist so all args are available.
      RemdOutArgs->Reset();
      if ( REMDtrajout.Add(NULL, RemdOutArgs, trajParm) ) {
        mprinterr("    Error: remdtraj remdout: Could not set up output traj %s\n",repFilename);
        j=-1;
      }
      // Close input traj
      (*it)->EndTraj();
      // If j==-1 an error occured, bail out
      if (j==-1) break; 
    } // END LOOP over input remd trajectories
    free(repFilename);
    free(Prefix);
    // If j==-1 an error occurred.
    if (j==-1) return 1;
  } // END REMDtrajout

  return 0;
}

/* RemdTraj::BeginTraj()
 * Open each trajectory in the list. 
 */
int RemdTraj::BeginTraj(bool showProgress) {
  std::list<TrajectoryFile *>::iterator traj;
  // Open the trajectory
  for (traj=REMDtraj.begin(); traj!=REMDtraj.end(); traj++)
    (*traj)->BeginTraj(false);
  numFramesRead=0;

  // Set up a progress bar
  if (showProgress) progress = new ProgressBar(stop);

  rprintf( "-REMD [%s] (%i-%i, %i) -----\n",trajName,start+1,stop+1,offset);

  // Open output trajectories
  if (RemdOutArgs!=NULL) {
    for (traj=REMDtrajout.begin(); traj!=REMDtrajout.end(); traj++)
      (*traj)->BeginTraj(false);
  }

  return 0;
}

/* RemdTraj::EndTraj()
 * Close all trajectories in the list.
 */
int RemdTraj::EndTraj() {
  std::list<TrajectoryFile *>::iterator traj;

  for (traj=REMDtraj.begin(); traj!=REMDtraj.end(); traj++)
    (*traj)->EndTraj();
  // Close output trajectories
  if (RemdOutArgs!=NULL) {
    for (traj=REMDtrajout.begin(); traj!=REMDtrajout.end(); traj++)
      (*traj)->EndTraj();
  }
  return 0;
}

/* RemdTraj::GetNextFrame()
 */
int RemdTraj::GetNextFrame(double *X, double *box, double *T) {
  std::list<TrajectoryFile *>::iterator reptrajin;
  std::list<TrajectoryFile *>::iterator reptrajout;
  int trajout, nrep;
  bool found = false;
  
  for (reptrajin=REMDtraj.begin(); reptrajin!=REMDtraj.end(); reptrajin++) {
    // No conversion to replica trajectories: Just find target temp
    if (RemdOutArgs==NULL) {
      if ((*reptrajin)->GetNextFrame(X,box,T)) return 1;
      // Check if this is the target temp
      if (*T == remdtrajtemp) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        return 0;
      }

    // All input REMD trajectories converted to temperature trajectories. 
    } else {
      // remdframe is allocated in SetupRead
      if ((*reptrajin)->GetNextFrame(remdframe->X,remdframe->box,&(remdframe->T))) return 1;
      // Check if this is the target temp. If so, set main Frame coords/box/temp
      if (remdframe->T == remdtrajtemp) {
        //printf("REMDTRAJ: remdout: Set %i TEMP=%lf\n",set,remdframe->T);
        for (int x=0; x < remdframe->N; x++)
          X[x] = remdframe->X[x];
        for (int b=0; b < 6; b++)
          box[b] = remdframe->box[b];
        *T = remdframe->T;
        found=true;
      }
      // Figure out which output file matches this temperature
      trajout=-1;
      nrep=0;
      for (reptrajout=REMDtrajout.begin(); reptrajout!=REMDtrajout.end(); reptrajout++) {
        if (TemperatureList[nrep] == remdframe->T) {
          trajout = nrep;
          break;
        }
        nrep++;
      }
      if (trajout==-1) {
        mprinterr("\nError: RemdTraj: remdout: Temperature %6.2lf not found in Temperature list.\n",
                  remdframe->T);
        return 1;
      }
      // Write to output traj
      (*reptrajout)->WriteFrame((*reptrajin)->CurrentFrame(),trajParm,
                                remdframe->X,remdframe->box,remdframe->T);
    }
  }  // END LOOP over input remd trajectories
  if (found) return 0;
  // If we have made it here this means target was not found
  mprinterr("\nError: RemdTraj: Final repTemp value read= %lf, set %i\n",*T,
            REMDtraj.front()->CurrentFrame());
  mprinterr("Could not find target %lf in any of the replica trajectories.\n",
          remdtrajtemp);
  mprinterr("Check that all replica trajectory files were found and that\n");
  mprinterr("none of the trajectories are corrupted (e.g. missing a temperature).\n");
  return 1;
}

/* RemdTraj::PrintInfo()
 */
void RemdTraj::PrintInfo(int showExtended) {
  mprintf("REMD trajectories (%i total, lowest replica: %s)\n",numReplicas,trajName);
  REMDtraj.Info(8);
  if (RemdOutArgs!=NULL) {
    mprintf("        remdout: trajectories will be converted to temperature trajectories:\n");
    REMDtrajout.Info(8);
  }
  mprintf("        Looking for frames at %8.2lf K",remdtrajtemp);
}
    
