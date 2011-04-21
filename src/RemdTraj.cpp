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
}

// DESTRUCTOR
RemdTraj::~RemdTraj() {
  //fprintf(stderr,"RemdTraj Destructor.\n");
  if (replicaName!=NULL) free(replicaName);
}

/*
 * RemdTraj::SetReplicaName()
 * Set lowest replica name. Also set target temperature.
 */
void RemdTraj::SetReplicaName(char* repName, ArgList *A) {
  replicaName=(char*) malloc( (strlen(repName)+1) * sizeof(char) );
  strcpy(replicaName,repName); 
  remdtrajtemp=A->getKeyDouble("remdtrajtemp",0.0);
}

/*
 * RemdTraj::NoTempInfo()
 * Check that trajectory has temperature info.
 * Return 1 if no temperature info, 0 otherwise.
 */
int RemdTraj::NoTempInfo(TrajFile *T) {

  if (T==NULL) return 1;

  if (T->hasTemperature==0) {
    mprintf("Error: REMDTRAJ: Trajectory %s does not contain temperature information.\n",
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
    mprintf("    ERROR: REMDTRAJ: Could not set up lowest replica file %s\n",replicaName);
    return 1;
  }
 
  // Set up lowest replica trajectory
  T->P = P;
  if ( T->SetupRead() ) {
    mprintf("    ERROR: REMDTRAJ: Setting up lowest replica %s for read.\n",replicaName);
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

  /*
   *  Scan for additional REMD traj files.
   *  Assume the extension of given trajectory is the number
   *  of the lowest replica, and that the other files are in
   *  sequence (e.g. rem.000, rem.001, rem.002 etc).
   */
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
    mprintf("    ERROR: REMDTRAJ: Could not find numeric extension.\n");
    mprintf("           Check that REMD files have naming scheme NAME.X\n");
    mprintf("           where X is an integer of arbitrary width.\n");
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
      mprintf("    ERROR: REMDTRAJ: Character #%i (%c) in extension %s is not a number!\n",
              j,ReplicaExt[j],ReplicaExt);
      free(ReplicaExt);
      free(Prefix);
      return 1;
    }
  }

  // Store lowest replica number
  j=atoi(ReplicaExt);
  //mprintf("  REMDDEBUG: index of first replica = %i\n",j);

  // Assume replica file names all have same length
  repFilename=(char*) malloc( (strlen(T->File->basefilename)+1) * sizeof(char));

  /*
   * Search for a replica number lower than this. Correct functioning
   * of the replica code requires the file specified by trajin be the
   * lowest # replica.
   */
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
      mprintf("    ERROR: REMDTRAJ: Setting up replica %s for read.\n",repFilename);
      delete T;
      free(repFilename);
      free(Prefix);
      free(ReplicaExt);
      return 1;
    }
    // Check for temperature information
    if ( NoTempInfo(T) ) {
      delete T;
      free(repFilename);
      free(Prefix);
      free(ReplicaExt);
      return 1;
    }
    // Check that #Frames and box info matches
    if ( Frames!=T->Frames || BoxType!=T->BoxType ) {
      mprintf(
              "    ERROR: REMDTRAJ: #Frames (%i) or box type (%i) in replica does not match\n",
              T->Frames, T->BoxType);
      mprintf("                     values in lowest replica (Frames=%i, boxtype=%i)\n",
              Frames,BoxType);
      delete T;
      free(repFilename);
      free(Prefix);
      free(ReplicaExt);
      return 1;
    }
    // Add it to the list
    REMDtraj.push_back(T);
    // Increment
    j++;
    numReplicas++;
    sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,j,CompressExt);
  }
  mprintf("REMDTRAJ:\n");
  REMDtraj.Info();
  
  free(repFilename);
  free(Prefix);
  free(ReplicaExt);
  return 0;
}

int RemdTraj::open() {
//  TrajFile *T;
  std::list<TrajFile *>::iterator it;

  for (it=REMDtraj.begin(); it!=REMDtraj.end(); it++)
    (*it)->Begin();

//  REMDtraj.Begin();
//  while ( (T=REMDtraj.Next())!=NULL ) {
//    // No allocation etc, just open
//    T->Begin(); 
//  }
  return 0;
}

void RemdTraj::close() {
  std::list<TrajFile *>::iterator it;

  for (it=REMDtraj.begin(); it!=REMDtraj.end(); it++)
    (*it)->End();

//  TrajFile *T;
  
//  REMDtraj.Begin();
//  while ( (T=REMDtraj.Next())!=NULL ) {
//    T->End();
//  }
}

/*
 * RemdTraj::getFrame()
 */
int RemdTraj::getFrame(int set) {
//  TrajFile *T;
  std::list<TrajFile *>::iterator it;
  
  for (it=REMDtraj.begin(); it!=REMDtraj.end(); it++) {
//  REMDtraj.Begin();
//  while ( (T=REMDtraj.Next())!=NULL ) {
    // Set replica frame to this frame (F should have been allocated already)
    (*it)->F=F;
    if ((*it)->getFrame(set)) return 1;
    // F belongs to remdtraj, not the individual replica trajectory so no need to free
    (*it)->F=NULL;
    // Check if this is the target temp
    if (F->T == remdtrajtemp) {
      //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
      return 0;
    }
  }
  // If we have made it here this means target was not found
  mprintf("\nREMDTRAJ: Final repTemp value read= %lf, set %i\n",F->T,set);
  mprintf("Could not find target %lf in any of the replica trajectories.\n",
          remdtrajtemp);
  mprintf("Check that all replica trajectory files were found and that\n");
  mprintf("none of the trajectories are corrupted (e.g. missing a temperature).\n");
  return 1;
}

/*
 * Info()
 */
void RemdTraj::Info() {
  mprintf("REMD trajectories (%i total, lowest replica: %s)\n",numReplicas,replicaName);
  mprintf("        Looking for frames at %8.2lf K",remdtrajtemp);
}
    
