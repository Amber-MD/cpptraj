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
  repOutname=NULL;
  lowestRepnum=0;

  hasTrajout=false;
  remdtrajtemp=0.0;
  TemperatureList=NULL;
  // Used for writing replica trajectories
  remdX = NULL;
  remdV = NULL;
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
  if (repOutname!=NULL) free(repOutname);
  if (CompressExt!=NULL) free(CompressExt);
  if (remdX!=NULL) free(remdX);
  if (remdV!=NULL) free(remdV);
}

// RemdTraj::SetupTemperatureList()
/** For each REMD input traj get the temperature of the first frame to
  * create a complete list of replica temperatures. This list is used
  * when writing out temperature trajectories from input replica
  * trajectories.
  */
int RemdTraj::SetupTemperatureList(int natom) {
  std::vector<TrajectoryIO*>::iterator replica;
  int repnum = 0;
  int err = 0;

  // Allocate space for reading in coords and velocities
  remdN = natom * 3; 
  remdX = (double*) malloc( remdN * sizeof(double));
  if (remdX==NULL) return 1;
  if (hasVelocity) {
    remdV = (double*) malloc( remdN * sizeof(double));
    if (remdV==NULL) return 1;
  }
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
      err = ((*replica)->readFrame(0,remdX, remdV, remdbox, TemperatureList+repnum));
      mprintf("      Rep %i T = %6.2lf\n",repnum,TemperatureList[repnum]);
      repnum++;
      // Close input traj
      (*replica)->closeTraj();
    }
    // If err!=0 an error occurred, bail out.
    if (err!=0) break;
  }
  if (err==0) hasTrajout=true;
  return err;
}

// RemdTraj::SetReplicaName()
/** Given the name of the lowest replica file, set basic replica name 
  * information. Get the file prefix up to the numerical extension, the 
  * width of the numerical extension, and the compression extension if 
  * the file is compressed.
  * Set and return the number of the lowest replica, or -1 on error.
  */
int RemdTraj::SetReplicaName(char* filenameIn) {
  CpptrajFile remdfile;
  int lastChar,lastDot,j;
  char *ReplicaExt;

  // STEP 1 - Get filename Prefix
  // Assume the extension of given trajectory is the number of the lowest 
  // replica, and that the other files are in sequence (e.g. rem.000, rem.001, 
  // rem.002 or rem.000.gz, rem.001.gz, rem.002.gz etc).
  // Set up a CpptrajFile here to check compression and extension.  
  if (remdfile.SetupFile(filenameIn, READ, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug)) 
    return -1;
  if (remdfile.compressType!=NO_COMPRESSION && remdfile.Ext!=NULL) {
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

// RemdTraj::GetReplicaName()
/** If name information has already been set by SetRepicaName, return
  * the expected replica filename given a replica number.
  */
char *RemdTraj::GetReplicaName(int repnum) {
  if (repFilename==NULL) return NULL;
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,repnum,CompressExt);
  return repFilename;
}

// RemdTraj::GetLowestReplicaName()
/** If name information has already been set by SetRepicaName, return
  * the lowest replica filename.
  */
// NOTE: obsolete?
char *RemdTraj::GetLowestReplicaName() {
  if (repFilename==NULL) return NULL;
  sprintf(repFilename,"%s.%0*i%s",Prefix,ExtWidth,lowestRepnum,CompressExt);
  return repFilename;
}

// RemdTraj::GetTemperatureName()
/** If the temperature list has been set up, return a filename with the 
  * specified replica numbers temperature appended to it. Note that the
  * replica numbers in this case start from 0 as opposed to incoming
  * replica filenames, where the lowest replica number can be arbitrary.
  */
char *RemdTraj::GetTemperatureName(char *prefix, int repnum) {
  if (TemperatureList==NULL) {
    mprinterr("Error: RemdTraj::GetTemperatureName: TemperatureList not set.\n");
    return NULL;
  }
  if (prefix==NULL) return NULL;
  if (repnum<0 || repnum > (int) REMDtraj.size()) return NULL;
  repOutname = (char*) realloc(repOutname, (strlen(prefix)+32) * sizeof(char));
  sprintf(repOutname,"%s.%6.2lf",prefix,TemperatureList[repnum]);
  return repOutname;
}

// RemdTraj::openTraj()
/** Open each trajectory in the list. */
int RemdTraj::openTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // DEBUG
  mprintf("REMD: OPENING REMD TRAJECTORIES\n");
  // Open the trajectory
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++)
    if ((*replica)->openTraj()) return 1;

  // Open output trajectories for writing
  if (hasTrajout) {
    for (replica=REMDtrajout.begin(); replica!=REMDtrajout.end(); replica++) 
      if ((*replica)->openTraj()) return 1;
  }

  return 0;
}

// RemdTraj::closeTraj()
/** Close all trajectories in the list. */
void RemdTraj::closeTraj() {
  std::vector<TrajectoryIO *>::iterator replica;
  // Close the trajectory
  for (replica=REMDtraj.begin(); replica!=REMDtraj.end(); replica++)
    (*replica)->closeTraj();

  // Close output trajectories
  if (hasTrajout) {
    for (replica=REMDtrajout.begin(); replica!=REMDtrajout.end(); replica++)
      (*replica)->closeTraj();
  }
}

// RemdTraj::readFrame()
/** Read the next frame from the list of input trajectories. Choose the
  * one that matches remdtrajtemp. Write trajectories if specified.
  */
int RemdTraj::readFrame(int set, double *X, double *V, double *box, double *T) {
  std::vector<TrajectoryIO *>::iterator reptrajin;
  std::vector<TrajectoryIO *>::iterator reptrajout;
  int trajout, nrep;
  bool found = false;
  
  for (reptrajin=REMDtraj.begin(); reptrajin!=REMDtraj.end(); reptrajin++) {
    // No conversion to replica trajectories: Just find target temp
    if (!hasTrajout) {
      if ((*reptrajin)->readFrame(set,X,V,box,T)) return 1;
      // Check if this is the target temp
      if (*T == remdtrajtemp) {
        //printf("REMDTRAJ: Set %i TEMP=%lf\n",set,F->T);
        return 0;
      }

    // All input REMD trajectories converted to temperature trajectories. 
    } else {
      // remdX, remdV, remdN is allocated in SetupTemperatureList 
      if ((*reptrajin)->readFrame(set,remdX,remdV,remdbox,&remdT)) return 1;
      // Check if this is the target temp. If so, set main Frame coords/box/temp
      if (remdT == remdtrajtemp) {
        //mprintf("REMDTRAJ: remdout: Set %i TEMP=%lf\n",set+1,remdT);
        memcpy(X, remdX, remdN * sizeof(double));
        if (V!=NULL && hasVelocity) 
          memcpy(V, remdV, remdN*sizeof(double));
        memcpy(box, remdbox, 6 * sizeof(double));
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
      //mprintf("REMDTRAJ: Write set %i, T %6.2lf, to %s\n",set+1, remdT, 
      //        (*reptrajout)->Filename());
      (*reptrajout)->writeFrame(set, remdX,remdV,remdbox,remdT);
    }
  }  // END LOOP over input remd trajectories
  if (found) return 0;
  // If we have made it here this means target was not found
  mprinterr("\nError: RemdTraj: Final repTemp value read= %lf, set %i\n",*T,set+OUTPUTFRAMESHIFT);
  mprinterr("Could not find target %lf in any of the replica trajectories.\n",
            remdtrajtemp);
  mprinterr("Check that all replica trajectory files were found and that\n");
  mprinterr("none of the trajectories are corrupted (e.g. missing a temperature).\n");
  return 1;
}

// RemdTraj::info()
void RemdTraj::info() {
  mprintf("REMD trajectories (%i total)\n", (int)REMDtraj.size());
  if (hasTrajout) {
    mprintf("  remdout: trajectories will be converted to temperature trajectories (%i total):\n",
            (int)REMDtrajout.size());
    for (int repnum=0; repnum < (int)REMDtraj.size(); repnum++)
      mprintf("\t\t[%6.2lf]\n",TemperatureList[repnum]);
  }
  mprintf("        Looking for frames at %8.2lf K",remdtrajtemp);
}
 
