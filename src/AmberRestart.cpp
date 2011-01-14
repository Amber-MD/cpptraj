// AmberRestart
#include <cstdio> // sscanf, sprintf
#include <cstdlib>
#include <cstring>
#include "AmberRestart.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AmberRestart::AmberRestart() {
  frameSize=0;
  frameBuffer=NULL;
  numBoxCoords=0;
  hasVelocity=0;
  restartAtoms=0;
  restartTime=0;
  restartTemp=-1.0;
}

// DESTRUCTOR
AmberRestart::~AmberRestart() {
  //fprintf(stderr,"AmberRestart Destructor.\n");
  // Close file. Should be safe
  close();
  if (frameBuffer!=NULL) free(frameBuffer);
}

/*
 * AmberRestart::close()
 */
void AmberRestart::close() {
  File->CloseFile();
}

/*
 * AmberRestart::open()
 */
int AmberRestart::open() {
  int nread,titleSize;
  char buffer[83];

  switch (File->access) {

    case READ :
      if (File->OpenFile()) return 1;
      // Read in title
      if (title==NULL) title=(char*) malloc(83*sizeof(char));
      if (File->IO->Gets(title,82)) {
         mprintf("Error: AmberRestart::open(): Reading restart title.\n");
        return 1;
      }
      // Read in natoms, time, and Replica Temp if present
      if (File->IO->Gets(buffer,82)) {
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
      } else if (nread==2) {
        restartTemp=-1.0;
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

/*
 * AmberRestart::SetupWrite()
 */
int AmberRestart::SetupWrite() {
  int frame_lines;

  restartAtoms=P->natom;
  // Calculate the length of coordinate frame in bytes
  frame_lines = (restartAtoms * 3) / 6;
  if (((restartAtoms * 3) % 6) > 0)
    frame_lines++;
  frameSize = ((restartAtoms * 3 * 12) + frame_lines);

  // Dont know ahead of time if velocities will be used, allocate space
  // just in case
  frameSize+=frameSize;

  // If box coords are present, allocate extra space for them
  if (isBox>0) {
    numBoxCoords=6;
    frameSize+=((numBoxCoords*12)+1);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer=(char*) calloc(frameSize+1, sizeof(char));

  return 0;
}

/*
 * AmberRestart::SetupRead()
 * Set up amber restart file for reading. Check that number of atoms matches
 * number of atoms in associated parmtop. Check for box/velocity info.
 */
int AmberRestart::SetupRead() {
  char buffer[83];
  int frame_lines,lineSize;

  open(); // Gets title, time, natoms, and temp if present

  // Check that natoms matches parm natoms
  if (restartAtoms!=P->natom) {
    mprintf("Error: AmberRestart::SetupRead(): Number of atoms in restart (%i)  \n",
            restartAtoms);
    mprintf("       does not match number of atoms in parmtop (%i)\n",P->natom);
    return 1;
  }

  // Calculate the length of coordinate frame in bytes
  frame_lines = (restartAtoms * 3) / 6;
  if (((restartAtoms * 3) % 6) > 0)
    frame_lines++;
  // For DOS files CR present before newline
  if (File->isDos) frame_lines*=2;
  frameSize = ((restartAtoms * 3 * 12) + frame_lines);
  frameBuffer=(char*) malloc(frameSize*sizeof(char));
  //if (debug>0) mprintf("    Amber Restart frameSize= %i\n",frameSize);

  // Read past restart coords 
  if ( File->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    mprintf("Error: AmberRestart::SetupRead(): Error reading coordinates.\n");
    return 1; 
  }

  // Attempt a second read to get velocities or box coords
  lineSize = File->IO->Read(frameBuffer,sizeof(char),frameSize);

  // If 0 or -1 no box or velo 
//  if (lineSize==-1) {
  if (lineSize<=0) {
//    mprintf("Error: AmberRestart::SetupRead(): Error checking for box/velocity info.\n");
//    return 1;

  // If 0 probably at EOF. No box or velo.
//  } else if (lineSize==0) {
    isBox=0;
    hasVelocity=0;

  // If 36 or 72 (+1 newline) box info.
  } else if (lineSize==37 || lineSize==73) {
    isBox=1;
    numBoxCoords = (lineSize-1) / 12;
    hasVelocity=0;

  // If filled framebuffer again, has velocity info. Check for box after velocity.
  } else if (lineSize==frameSize) {
    hasVelocity=1;
    if (File->IO->Gets(buffer,82)==0) {
      lineSize=strlen(buffer);
      if (lineSize==37 || lineSize==73) {
        isBox=1;
        numBoxCoords = (lineSize-1) / 12;
      } else {
        mprintf("Error: AmberRestart::SetupRead():\n");
        mprintf("       Expect only 3 or 6 box coords in box coord line.\n");
        return 1;
      }
    } else
      isBox=0;

  // Otherwise, who knows what was read?
  } else {
    mprintf("Error: AmberRestart::SetupRead(): When attempting to read in\n");
    mprintf("       box coords/velocity info got %i chars, expected 0, 37,\n",lineSize);
    mprintf("       73, or %i.\n",frameSize);
    return 1;
  }

  // Recalculate the frame size
  if (hasVelocity) frameSize+=frameSize;
  if (isBox) frameSize+=( (numBoxCoords*12) + 1 );
  frameBuffer=(char*) realloc(frameBuffer, frameSize*sizeof(char));

  if (debug > 0) {
    mprintf("    Amber Restart isBox=%i hasVelocity=%i numBoxCoords=%i\n",
            isBox,hasVelocity,numBoxCoords);
    mprintf("    Amber Restart frameSize= %i\n",frameSize);
  }
  
  // Set TrajFile variables
  seekable=1; // Set seekable since only 1 frame read (i.e. we know size)
  Frames=1;
  stop=1;

  close();
  return 0;
}

/*
 * AmberRestart::getFrame()
 * Get the restart file frame. If velocities are present, read those too.
 */
int AmberRestart::getFrame(int set) {
  char *bufferPosition;
  // Read restart coords into frameBuffer
  if ( File->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    mprintf("Error: AmberRestart::getFrame(): Error reading coordinates.\n");
    return 1;
  }

  // Set frame temp
  F->T = restartTemp;

  // Convert coords to Frame - Frame expects N coords to be present in buffer
  if ( (bufferPosition = F->BufferToFrame(frameBuffer, 12))==NULL ) {
    mprintf("Error: AmberRestart::getFrame: * detected in coordinates of %s\n",trajfilename);
    return 1;  
  }
  // Convert velocity to Frame if present
  if (hasVelocity) {
    if (F->V==NULL) F->V = new Frame(restartAtoms,NULL);
    if ( (bufferPosition = F->V->BufferToFrame(bufferPosition, 12))==NULL ) {
      mprintf("Error: AmberRestart::getFrame: * detected in velocities of %s\n",
              trajfilename);
      return 1;
    }
    //F->V->printAtomCoord(0);
  }
  // Convert box to Frame if present
  if (isBox) {
    if ( (bufferPosition = F->BufferToBox(bufferPosition, numBoxCoords, 12))==NULL ) {
      mprintf("Error: AmberRestart::getFrame: * detected in box coordinates of %s\n",
              trajfilename);
      return 1;
    }
  }

  return 0;
}

/*
 * AmberRestart::writeFrame()
 * Write coords in Frame to file in amber restart format. Always calculate the 
 * frame size since coords may have been stripped from Frame.
 */
int AmberRestart::writeFrame(int set) {
  char buffer[1024];
  char *bufferPosition;

  sprintf(buffer,"%s.%i",File->filename,set);
  if (File->IO->Open(buffer,"wb")) return 1;

  // Write out title
  File->IO->Printf("%-80s\n",title);
  // Write out atoms, time, and temp (if not -1)
  // NOTE: Use F->natom instead of restart Atoms in case of strip cmd
  //restartTime = (double) set; // NOTE: Could be an option to use set eventually
  File->IO->Printf("%5i%15.7lE\n",F->natom,restartTime);

  // Write coords to buffer
  bufferPosition = F->FrameToBuffer(frameBuffer,"%12.7lf",12,6);
  // Write velocity to buffer
  if (F->V!=NULL)  // NOTE: Use hasVelocity in addition/instead?
    bufferPosition = F->V->FrameToBuffer(bufferPosition,"%12.7lf",12,6);
  // Write box to buffer
  if (isBox)
    bufferPosition = F->BoxToBuffer(bufferPosition, numBoxCoords, "%12.7lf",12);

  //if (seekable) fseek(fp, titleSize+(set*frameSize),SEEK_SET);
  frameSize = (int) (bufferPosition - frameBuffer);

  if (File->IO->Write(frameBuffer,sizeof(char),frameSize)) return 1;

  File->IO->Close();

  return 0;
}

/*
 * Info()
 */
void AmberRestart::Info() {
  mprintf("  File (%s) is an AMBER restart file", File->filename);
  if (hasVelocity) mprintf(" with velocity info");
}
