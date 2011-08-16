// AmberRestart
#include <cstdio> // sscanf
#include <cstdlib>
#include <cstring>
#include "Traj_AmberRestart.h"
#include "CpptrajStdio.h"
#include "CharBuffer.h"

// CONSTRUCTOR
AmberRestart::AmberRestart() {
  restartAtoms=0;
  frameSize=0;
  frameBuffer=NULL;
  numBoxCoords=0;
  restartTime=0;
  restartTemp=-1.0;
  // Set seekable since only 1 frame read (i.e. we know size) 
  seekable=true;
}

// DESTRUCTOR
AmberRestart::~AmberRestart() {
  //fprintf(stderr,"AmberRestart Destructor.\n");
  // Close file. Should be safe
  closeTraj();
  if (frameBuffer!=NULL) free(frameBuffer);
}

/* AmberRestart::closeTraj()
 */
void AmberRestart::closeTraj() {
  tfile->CloseFile();
}

/* AmberRestart::openTraj()
 * Open the restart file. Get title, time, restart atoms, temperature
 */
int AmberRestart::openTraj() {
  int nread,titleSize;
  char buffer[83];

  switch (tfile->access) {

    case READ :
      if (tfile->OpenFile()) return 1;
      // Read in title
      if (tfile->IO->Gets(buffer,82)) {
         mprintf("Error: AmberRestart::open(): Reading restart title.\n");
        return 1;
      }
      SetTitle(buffer);
      // Read in natoms, time, and Replica Temp if present
      if (tfile->IO->Gets(buffer,82)) {
        mprintf("Error: AmberRestart::open(): Reading restart atoms/time.\n");
        return 1;
      }
      nread=sscanf(buffer,"%5i%15lE%15lE",&restartAtoms,&restartTime,&restartTemp);
      if (nread<1) {
        mprintf("Error: AmberRestart::open(): Getting restart atoms/time.\n");
        return 1;
      } else if (nread==1) {
        restartTime=0.0;
        restartTemp=-1.0;
        hasTemperature=false;
      } else if (nread==2) {
        restartTemp=-1.0;
        hasTemperature=false;
      } else {
        hasTemperature=true;
      }
      if (debug>0) mprintf("  Amber restart: Atoms=%i Time=%lf Temp=%lf\n",restartAtoms,
                           restartTime, restartTemp);
      break;

    case APPEND :
      mprintf("Error: Append not supported for Amber Restart files.\n");
      return 1;
      break;

    case WRITE :
      // Set up title
      if (title==NULL) {
        title=(char*) malloc(81*sizeof(char));
        strcpy(title,"Cpptraj Generated Restart");
        titleSize=81;
      } else {
        titleSize=strlen(title);
        if (titleSize>81) {
          mprintf("Error: AmberTraj::open: Given title is too long!\n[%s]\n",title);
          return 1;
        }
      }

  }

  return 0; 
}

/* AmberRestart::SetNoVelocity()
 */
void AmberRestart::SetNoVelocity() {
  hasVelocity=false;
}

/* AmberRestart::processWriteArgs()
 */
int AmberRestart::processWriteArgs(ArgList *argIn) {
  // For write, assume we want velocities unless specified
  hasVelocity=true;
  if (argIn->hasKey("novelocity")) this->SetNoVelocity();
  return 0;
}

/* AmberRestart::setupWrite()
 * Allocate a character buffer based on number of coords and whether 
 * velocities/box info is present.
 */
int AmberRestart::setupWrite(AmberParm *trajParm) {
  int frame_lines;

  restartAtoms=trajParm->natom;
  natom3=trajParm->natom * 3;
  // Calculate the length of coordinate frame in bytes
  frame_lines = (natom3) / 6;
  if ((natom3 % 6) > 0)
    frame_lines++;
  frameSize = ((natom3 * 12) + frame_lines);

  // Dont know ahead of time if velocities will be used, allocate space
  // just in case. Veolcity will not be written if V input is NULL.
  frameSize+=frameSize;

  // If box coords are present, allocate extra space for them
  if (hasBox) {
    numBoxCoords=6;
    frameSize+=((numBoxCoords*12)+1);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer=(char*) calloc(frameSize+1, sizeof(char));

  return 0;
}

/* AmberRestart::getBoxAngles()
 * Based on input buffer, determine num box coords and get box angles.
 */
void AmberRestart::getBoxAngles(char *boxline, int boxlineSize) {
  double box[6];
  numBoxCoords = (boxlineSize-1) / 12;
  if (numBoxCoords>3) {
    sscanf(boxline, "%8lf%8lf%8lf%8lf%8lf%8lf",box,box+1,box+2,boxAngle,boxAngle+1,boxAngle+2);
  } 
}

/* AmberRestart::setupRead()
 * Set up amber restart file for reading. Check that number of atoms matches
 * number of atoms in associated parmtop. Check for box/velocity info.
 */
int AmberRestart::setupRead(AmberParm *trajParm) {
  char buffer[83];
  int frame_lines,lineSize;

  if (openTraj()) return -1; // Gets title, time, natoms, and temp if present

  // Check that natoms matches parm natoms
  if (restartAtoms!=trajParm->natom) {
    mprinterr("Error: Number of atoms in Amber Restart %s (%i) does not\n",
              tfile->filename, restartAtoms);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->natom);
    return -1;
  }
  natom3 = restartAtoms * 3;

  // Calculate the length of coordinate frame in bytes
  frame_lines = natom3 / 6;
  if ((natom3 % 6) > 0)
    frame_lines++;
  // For DOS files CR present before newline
  if (tfile->isDos) frame_lines*=2;
  frameSize = ((natom3 * 12) + frame_lines);
  frameBuffer=(char*) malloc(frameSize*sizeof(char));
  //if (debug>0) mprintf("    Amber Restart frameSize= %i\n",frameSize);

  // Read past restart coords 
  if ( tfile->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    mprintf("Error: AmberRestart::SetupRead(): Error reading coordinates.\n");
    return -1; 
  }

  // Attempt a second read to get velocities or box coords
  lineSize = tfile->IO->Read(frameBuffer,sizeof(char),frameSize);

  // If 0 or -1 no box or velo 
  if (lineSize<=0) {
    hasBox=false; 
    hasVelocity=false;

  // If 36 or 72 (+1 newline) box info.
  } else if (lineSize==37 || lineSize==73) {
    hasBox=true;
    getBoxAngles(frameBuffer,lineSize);
    hasVelocity=false;

  // If filled framebuffer again, has velocity info. Check for box after velocity.
  } else if (lineSize==frameSize) {
    hasVelocity=true;
    if (tfile->IO->Gets(buffer,82)==0) {
      lineSize=strlen(buffer);
      if (lineSize==37 || lineSize==73) {
        hasBox=true;
        getBoxAngles(buffer,lineSize);
      } else {
        mprintf("Error: AmberRestart::SetupRead():\n");
        mprintf("       Expect only 3 or 6 box coords in box coord line.\n");
        return -1;
      }
    } else
      hasBox=false;

  // Otherwise, who knows what was read?
  } else {
    mprintf("Error: AmberRestart::SetupRead(): When attempting to read in\n");
    mprintf("       box coords/velocity info got %i chars, expected 0, 37,\n",lineSize);
    mprintf("       73, or %i.\n",frameSize);
    return -1;
  }

  // Recalculate the frame size
  if (hasVelocity) frameSize+=frameSize;
  if (hasBox) frameSize+=( (numBoxCoords*12) + 1 );
  frameBuffer=(char*) realloc(frameBuffer, frameSize*sizeof(char));

  if (debug > 0) {
    mprintf("    Amber Restart hasBox=%i hasVelocity=%i numBoxCoords=%i\n",
            (int)hasBox,(int)hasVelocity,numBoxCoords);
    mprintf("    Amber Restart frameSize= %i\n",frameSize);
  }
 
  closeTraj();
  // Only 1 frame in restart by definition
  return 1;
}

/* AmberRestart::readFrame()
 * Get the restart file frame. If velocities are present, read those too.
 */
int AmberRestart::readFrame(int set,double *X,double *V,double *box, double *T) {
  char *bufferPosition;
  // Read restart coords into frameBuffer
  if ( tfile->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    mprintf("Error: AmberRestart::getFrame(): Error reading coordinates.\n");
    return 1;
  }

  // Set frame temp
  if (hasTemperature)
    *T = restartTemp;

  // Get coords from buffer
  if ( (bufferPosition = BufferToDouble(frameBuffer, X, natom3, 12))==NULL ) {
    mprintf("Error: AmberRestart::getFrame: * detected in coordinates of %s\n",tfile->filename);
    return 1;  
  }
  // Get velocity from buffer if present
  if (hasVelocity && V!=NULL) {
    if ( (bufferPosition = BufferToDouble(bufferPosition, V, natom3, 12))==NULL ) {
      mprintf("Error: AmberRestart::getFrame: * detected in velocities of %s\n",
              tfile->filename);
      return 1;
    }
    //F->V->printAtomCoord(0);
  }
  // Get box from buffer if present
  if (hasBox) {
    if ( (bufferPosition = BufferToDouble(bufferPosition, box, numBoxCoords, 12))==NULL ) {
      mprintf("Error: AmberRestart::getFrame: * detected in box coordinates of %s\n",
              tfile->filename);
      return 1;
    }
  }

  return 0;
}

/* AmberRestart::writeFrame()
 * Write coords in Frame to file in amber restart format. Always calculate the 
 * frame size since coords may have been stripped from Frame.
 */
int AmberRestart::writeFrame(int set, double *X, double *V, double *box, double T) {
  char buffer[1024];
  char *bufferPosition;

  NumberFilename(buffer,tfile->filename,set + OUTPUTFRAMESHIFT);
  if (tfile->IO->Open(buffer,"wb")) return 1;

  // Write out title
  tfile->IO->Printf("%-80s\n",title);
  // Write out atoms, time, and temp (if not -1)
  // NOTE: Use F->natom instead of restart Atoms in case of strip cmd
  restartTime = (double) set; // NOTE: Could be an option to use set eventually
  tfile->IO->Printf("%5i%15.7lE",restartAtoms,restartTime + OUTPUTFRAMESHIFT);
  if (hasTemperature)
    tfile->IO->Printf("%15.7lE",T);
  tfile->IO->Printf("\n");

  // Write coords to buffer
  bufferPosition = DoubleToBuffer(frameBuffer,X,natom3,"%12.7lf",12,6);
  // Write velocity to buffer. Check V since velocity not known ahead of time
  if (hasVelocity && V!=NULL)
    bufferPosition = DoubleToBuffer(bufferPosition,V,natom3,"%12.7lf",12,6);
  //if (F->V!=NULL)  // NOTE: Use hasVelocity in addition/instead?
  //  bufferPosition = DoubleToBuffer(bufferPosition,F->V->X,F->N,"%12.7lf",12,6);
  // Write box to buffer
  if (hasBox)
    bufferPosition = BoxToBuffer(bufferPosition, box, numBoxCoords, "%12.7lf",12);

  //if (seekable) fseek(fp, titleSize+(set*frameSize),SEEK_SET);
  frameSize = (int) (bufferPosition - frameBuffer);

  if (tfile->IO->Write(frameBuffer,sizeof(char),frameSize)) return 1;

  // Since when writing amber restarts a number is appended to the filename
  // dont use the CloseFile function in PtrajFile, just close.
  tfile->IO->Close();

  return 0;
}

/* AmberRestart::info()
 */
void AmberRestart::info() {
  mprintf("is an AMBER restart file");
  if (hasVelocity) mprintf(" with velocity info");
}
