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
  replicaName=NULL;
  numReplicas=0;
  RemdOutArgs=NULL;
  TemperatureList=NULL;
  remdframe=NULL;
}

// DESTRUCTOR
RemdTraj::~RemdTraj() {
  //fprintf(stderr,"RemdTraj Destructor.\n");
  if (replicaName!=NULL) free(replicaName);
  if (RemdOutArgs!=NULL) delete RemdOutArgs;
  if (TemperatureList!=NULL) free(TemperatureList);
  if (remdframe!=NULL) delete remdframe;
}

/* RemdTraj::SetReplicaName()
 * Set lowest replica name. Also set target temperature. If the remdout
 * keyword is specified, split the argument list at the keyword and treat
 * all arguments following remdout as trajout arguments. Temperature
 * trajectories will be written based on those arguments.
 */
void RemdTraj::SetReplicaName(char* repName, ArgList *A) {

  if (repName==NULL) {
    mprinterr("Error: RemdTraj::SetReplicaName: No filename specified.\n");
    return;
  }
  replicaName=(char*) malloc( (strlen(repName)+1) * sizeof(char) );
  strcpy(replicaName,repName); 
  remdtrajtemp=A->getKeyDouble("remdtrajtemp",0.0);
  if (A->hasKey("remdout")) {
    RemdOutArgs = A->SplitAt("remdout");
  }
}

/* RemdTraj::NoTempInfo()
 * Check that trajectory has temperature info.
 * Return 1 if no temperature info, 0 otherwise.
 */
int RemdTraj::NoTempInfo(TrajFile *T) {

  if (T==NULL) return 1;

  if (T->hasTemperature==0) {
    mprinterr("Error: REMDTRAJ: Trajectory %s does not contain temperature information.\n",
            T->File->filename);
    return 1;
  }

  return 0;
}

/*
 * RemdTraj::SetupRead()
 * Set up this trajectory for processing as an REMD trajectory. All other
 * replica trajectory files are searched for, and the base filename, starting
 * replica, and number of replicas is recorded so the files may be opened
 * later during processing.
 */
int RemdTraj::SetupRead() {
  int lastChar;     // The last character before compress extension (if any)
  int i,j, ExtWidth;
  TrajFile *T;
  char *ReplicaExt;  // Hold numeric extension
  char *CompressExt; // Hold compression extension
  char *Prefix;      // Hold filename up until the numeric extension
  char *repFilename; // Additional replica filenames

  if (debug>0) {
    mprintf("    REMDTRAJ: Using specified file as lowest replica: %s\n",replicaName);
    mprintf("    REMDTRAJ: Frames at %lf K will be processed.\n",remdtrajtemp);
  }

  // Set up lowest replica file
  T = REMDtraj.SetupTrajectory(replicaName,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE);
  if (T==NULL) {
    mprinterr("    ERROR: REMDTRAJ: Could not set up lowest replica file %s\n",replicaName);
    return 1;
  }
 
  // Set up lowest replica trajectory
  T->P = P;
  if ( T->SetupRead() ) {
    mprinterr("    ERROR: REMDTRAJ: Setting up lowest replica %s for read.\n",replicaName);
    delete T;
    return 1;
  }
  // Check for temperature information
  if ( NoTempInfo(T) ) return 1;
  // Use this trajectory to set up overall stop, frames, and box info.
  // NOTE: Should check that this is the case for ALL frames.
  stop = T->Frames;
  Frames = T->Frames;
  BoxType = T->BoxType;
  trajfilename = T->File->basefilename;
  // Add it to the list
  REMDtraj.push_back(T);
  //REMDtraj.Info();
  numReplicas=1;

  // Scan for additional REMD traj files.
  // Assume the extension of given trajectory is the number
  // of the lowest replica, and that the other files are in
  // sequence (e.g. rem.000, rem.001, rem.002 etc).
  if (T->File->compressType!=NONE && T->File->Ext!=NULL) {
    CompressExt=T->File->Ext; 
    lastChar = strlen(T->File->basefilename) - strlen(T->File->Ext);
  } else {
    CompressExt=(char*)"";
    lastChar = strlen(T->File->basefilename);
  }

  // Find location of last '.' (not including compression extension) and store it in i
  i=-1;
  for (j=0; j<lastChar; j++)
    if (T->File->basefilename[j]=='.') i=j;
  if (i==-1) {
    mprinterr("    ERROR: REMDTRAJ: Could not find numeric extension.\n");
    mprinterr("           Check that REMD files have naming scheme NAME.X\n");
    mprinterr("           where X is an integer of arbitrary width.\n");
    return 1;
  }

  // Store filename up until the numeric extension
  Prefix=(char*) malloc( (i+1) * sizeof(char));
  strncpy(Prefix,T->File->basefilename,i);
  Prefix[i]='\0';
  //mprintf("  REMDDEBUG: Replica filename prefix: %s\n",Prefix);

  ExtWidth=lastChar - i - 1;
  //mprintf("  REMDDEBUG: Last . in %s located at %i\n",T->File->basefilename,i);
  //mprintf("  REMDDEBUG: Allocating %i for extension\n",ExtWidth+1);
  //mprintf("  REMDDEBUG: EXTwidth=%i\n",ExtWidth);
  ReplicaExt=(char*) malloc( (ExtWidth+1) * sizeof(char));
  strncpy(ReplicaExt, T->File->basefilename + i + 1, ExtWidth);
  ReplicaExt[ExtWidth]='\0'; 
  //mprintf("  REMDDEBUG: Replica extension is %s\n",ReplicaExt);

  // Check that all digits in extension are numbers 
  for (j=0; j<ExtWidth; j++) {
    if (isdigit(ReplicaExt[j])==0) {
      mprinterr("    ERROR: REMDTRAJ: Character #%i (%c) in extension %s is not a number!\n",
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

  // Assume replica file names all have same length
  repFilename=(char*) malloc( (strlen(T->File->basefilename)+1) * sizeof(char));

  // Search for a replica number lower than this. Correct functioning
  // of the replica code requires the file specified by trajin be the
  // lowest # replica.
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j-1,CompressExt);
  T = REMDtraj.SetupTrajectory(repFilename,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE);
  if (T!=NULL) {
    mprintf(
            "    WARNING: REMDTRAJ: Replica# found lower than file specified with trajin!\n");
    mprintf("              (Found %s)\n",repFilename);
    mprintf("              trajin <file> remdtraj requires lowest # replica!\n");
    delete T;
  }

  // Search for and add all replicas higher than this.
  j++;
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j,CompressExt);
  while ( (T=REMDtraj.SetupTrajectory(repFilename,READ,UNKNOWN_FORMAT,UNKNOWN_TYPE))!=NULL ) {
    T->P = P;
    if ( T->SetupRead() ) {
      mprinterr("    ERROR: REMDTRAJ: Setting up replica %s for read.\n",repFilename);
      delete T;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Check for temperature information
    if ( NoTempInfo(T) ) {
      delete T;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Check that #Frames and box info matches
    if ( Frames!=T->Frames || BoxType!=T->BoxType ) {
      mprinterr(
              "    ERROR: REMDTRAJ: #Frames (%i) or box type (%i) in replica does not match\n",
              T->Frames, T->BoxType);
      mprinterr("                     values in lowest replica (Frames=%i, boxtype=%i)\n",
              Frames,BoxType);
      delete T;
      free(repFilename);
      free(Prefix);
      return 1;
    }
    // Add it to the list
    REMDtraj.push_back(T);
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
    remdframe = new Frame(P->natom, NULL);
    // Allocate space for temperature list
    TemperatureList = (double*) malloc(numReplicas*sizeof(double));
    j=0;
    // BEGIN LOOP over input remd trajectories
    // Get a list of all temperatures present in input REMDtraj
    for (std::list<TrajFile *>::iterator it=REMDtraj.begin(); it!=REMDtraj.end(); it++) {
      (*it)->Begin();
      // Read 1 frame to get temperature
      (*it)->F = remdframe; 
      (*it)->getFrame(0);
      TemperatureList[j++] = (*it)->F->T;
      if (debug>0) mprintf("    Rep %i T = %6.2lf\n",j-1,TemperatureList[j-1]);
      // Set up output filename for this temperature
      sprintf(repFilename,"%s.%6.2lf",Prefix,(*it)->F->T);
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
      if ( REMDtrajout.Add(RemdOutArgs, P) ) {
        mprinterr("    Error: remdtraj remdout: Could not set up output traj %s\n",repFilename);
        j=-1;
      }
      // Reset frame and close input traj
      (*it)->F = NULL;
      (*it)->End();
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

/* RemdTraj::open()
 * Open each trajectory in the list. Use Begin with no args so that
 * no allocation is performed.
 */
int RemdTraj::open() {
  std::list<TrajFile *>::iterator it;

  for (it=REMDtraj.begin(); it!=REMDtraj.end(); it++)
    (*it)->Begin();

  return 0;
}

/* RemdTraj::close()
 * Close all trajectories in the list.
 */
void RemdTraj::close() {
  std::list<TrajFile *>::iterator it;

  for (it=REMDtraj.begin(); it!=REMDtraj.end(); it++)
    (*it)->End();
}

/* RemdTraj::getFrame()
 */
int RemdTraj::getFrame(int set) {
  std::list<TrajFile *>::iterator reptrajin;
  std::list<TrajFile *>::iterator reptrajout;
  int trajout, nrep;
  bool found = false;
  
  for (reptrajin=REMDtraj.begin(); reptrajin!=REMDtraj.end(); reptrajin++) {
    // No conversion to replica trajectories: Just find target temp
    if (RemdOutArgs==NULL) {
      // Set replica frame to this frame (F should have been allocated already)
      (*reptrajin)->F=F;
      if ((*reptrajin)->getFrame(set)) return 1;
      // F belongs to remdtraj, not the individual replica trajectory so no need to free
      (*reptrajin)->F=NULL;
      // Check if this is the target temp
      if (F->T == remdtrajtemp) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        return 0;
      }

    // All input REMD trajectories converted to temperature trajectories. 
    } else {
      // Set replica frame to remdframe
      (*reptrajin)->F=remdframe;
      if ((*reptrajin)->getFrame(set)) return 1;
      (*reptrajin)->F=NULL;
      // Check if this is the target temp. If so, set main Frame coords/box/temp
      if (remdframe->T == remdtrajtemp) {
        //printf("REMDTRAJ: remdout: Set %i TEMP=%lf\n",set,remdframe->T);
        for (int x=0; x < remdframe->N; x++)
          F->X[x] = remdframe->X[x];
        for (int b=0; b < 6; b++)
          F->box[b] = remdframe->box[b];
        F->T = remdframe->T;
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
      // Perform output traj setup if needed
      if ((*reptrajout)->skip==0) {
        if (debug>0) rprintf("    Setting up %s for remd WRITE (%i atoms\n",
                           (*reptrajout)->trajfilename,(*reptrajout)->P->natom);
        if ((*reptrajout)->SetupWrite()) return 1;
        if ((*reptrajout)->Begin()) return 1;
        (*reptrajout)->skip=1;
      }
      // Set output traj frame to remdframe and write
      (*reptrajout)->F=remdframe;
      //fprintf(stdout,"DEBUG: %20s: Writing %i\n",(*reptrajout)->File->filename,set);
      if ((*reptrajout)->writeFrame(set)) return 1;
      // Dont want this F deallocd since it belongs to remdframe 
      (*reptrajout)->F=NULL;
    }
  }  // END LOOP over input remd trajectories
  if (found) return 0;
  // If we have made it here this means target was not found
  mprinterr("\nError: RemdTraj: Final repTemp value read= %lf, set %i\n",F->T,set);
  mprinterr("Could not find target %lf in any of the replica trajectories.\n",
          remdtrajtemp);
  mprinterr("Check that all replica trajectory files were found and that\n");
  mprinterr("none of the trajectories are corrupted (e.g. missing a temperature).\n");
  return 1;
}

/* Info()
 */
void RemdTraj::Info() {
  mprintf("REMD trajectories (%i total, lowest replica: %s)\n",numReplicas,replicaName);
  REMDtraj.Info(8);
  if (RemdOutArgs!=NULL) {
    mprintf("        remdout: trajectories will be converted to temperature trajectories:\n");
    REMDtrajout.Info(8);
  }
  mprintf("        Looking for frames at %8.2lf K",remdtrajtemp);
}
    
