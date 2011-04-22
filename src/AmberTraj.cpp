#include <cstdio> //sprintf
#include <cstdlib>
#include <cstring>
#include "AmberTraj.h"
#include "CpptrajStdio.h"
// AmberTraj

// Size of REMD header
#define REMD_HEADER_SIZE 42

// CONSTRUCTOR
AmberTraj::AmberTraj() {
  //fprintf(stderr,"AmberTraj Constructor.\n");
  titleSize=0;
  frameSize=0;
  frameBuffer=NULL;
  numBoxCoords=0;
  hasREMD=0;
}

// DESTRUCTOR
AmberTraj::~AmberTraj() {
  //fprintf(stderr,"AmberTraj Destructor.\n");
  if (frameBuffer!=NULL) free (frameBuffer);
}    

void AmberTraj::close() {
  File->CloseFile();
}

#define ROUTINE "AmberTraj::open()"
/* 
 * AmberTraj::open()
 * Open Amber Trajectory, read past title and set titleSize
 * Always read title so file pointer is always correctly positioned
 * to read a frame.
 */
int AmberTraj::open() {

  switch (File->access) {
    case READ:
      File->OpenFile();
      // Read in title, set size in bytes 
      if (title==NULL)  
        title=(char*) malloc(BUFFER_SIZE * sizeof(char));
      if ( File->IO->Gets(title,BUFFER_SIZE) ) {
        rprintf( "WARNING in %s: EOF encountered during reading of\n", ROUTINE);
        rprintf( "   title from (%s)\n", File->filename);
        return 1;
      }
      titleSize=strlen(title);
      break;

    case WRITE:
      File->OpenFile();
      // Write title
      // NOTE: THIS WILL SCREW UP IF OPEN CALLED MULTIPLE TIMES FOR WRITE!
      File->IO->Rank_printf(0,"%-*s\n",titleSize-1,title);
      break;

    case APPEND:
      // First, open the file with read access to set the title
      File->access=READ;
      if (open()) return 1;
      close();
      // Reset the access and open append
      File->access=APPEND;
      File->OpenFile();
      break;
  }

  return 0;
}
#undef ROUTINE

/*
 * AmberTraj::getFrame()
 * Read coordinates from Amber trajectory into Frame. 
 * NOTE: Precalculate the header, coord, and box offsets.
 */
int AmberTraj::getFrame(int set) {
  char Temp[9];
  char *bufferPosition;
  off_t offset;

  if (seekable) {
    offset = (off_t) set;
    offset *= (off_t) frameSize;
    offset += (off_t) titleSize;
    File->IO->Seek(offset);
  }

  if ( File->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) return 1;

  // Get REMD Temperature if present
  if (hasREMD>0) {
    strncpy(Temp,frameBuffer+33,8);
    Temp[8]='\0';
    F->T = atof(Temp);
    //if (debug>0) fprintf(stderr,"DEBUG: Replica T is %lf (%s)\n",F->T,Temp);
  }
  // Get Coordinates - hasREMD is size in bytes of REMD header
  if ( (bufferPosition = F->BufferToFrame(frameBuffer+hasREMD,8))==NULL ) {
    rprintf("Error: AmberTraj::getFrame: * detected in coordinates of %s\n",trajfilename);
    return 1;  
  } 
  if (BoxType!=0) { 
    if ( (bufferPosition = F->BufferToBox(bufferPosition,numBoxCoords,8))==NULL ) {
      rprintf("Error: AmberTraj::getFrame: * detected in box coordinates of %s\n",
              trajfilename);
      return 1;
    }
    // Set box angles to parmtop default if not read in
    if (numBoxCoords==3) {
      F->box[3] = P->Box[3];
      F->box[4] = P->Box[4];
      F->box[5] = P->Box[5];
    }
  }
  return 0;
}

/*
 * AmberTraj::writeFrame()
 * Write coordinates from Frame to frameBuffer. frameBuffer must be large
 * enough to accomodate all coords in F (handled by SetupWrite).
 * NOTE: The output frame size is calcd here - should it just be precalcd?
 */
int AmberTraj::writeFrame(int set) {
  int outFrameSize;
  char *bufferPosition;
  off_t offset;

  bufferPosition = frameBuffer;

  if (hasREMD) {
    sprintf(bufferPosition,"REMD  %8i %8i %8i %8.3lf\n",0,set+OUTPUTFRAMESHIFT,
            set+OUTPUTFRAMESHIFT,F->T);
    bufferPosition += hasREMD;
  }

  bufferPosition = F->FrameToBuffer(bufferPosition,"%8.3lf",8,10);
  if (BoxType!=0) 
    bufferPosition = F->BoxToBuffer(bufferPosition,numBoxCoords,"%8.3lf",8);

  outFrameSize = (int) (bufferPosition - frameBuffer);
  
  //if (seekable) 
  // NOTE: Seek only needs to happen when traj file changes
  //offset = (off_t) currentFrame;
  //offset *= (off_t) outFrameSize;
  //offset += (off_t) titleSize;
  //File->IO->Seek( offset);

  if (File->IO->Write(frameBuffer,sizeof(char),outFrameSize)) return 1;

  currentFrame++;
 
  return 0;
}

#define ROUTINE "AmberTraj::setup()"
/*
 * AmberTraj::SetupRead()
 * Setup opens the given file for read access, sets information. 
 */
int AmberTraj::SetupRead() {
  char buffer[BUFFER_SIZE];
  int frame_lines;
  int lineSize;
  long long int file_size, frame_size;
  double box[6]; // For checking box coordinates

  // Attempt to open the file. open() sets the title and titleSize
  if (open()) return 1;

  // Read 1 line to check for REMD header
  if ( File->IO->Gets(buffer, BUFFER_SIZE) ) {
    rprintf("Error: AmberTraj::SetupRead(): Reading second line of trajectory.\n");
    return 1;
  }
  if (strncmp(buffer,"REMD",4)==0 || strncmp(buffer,"HREMD",5)==0) {
    hasREMD=REMD_HEADER_SIZE;
    hasTemperature=1;
  }
  // Reopen the file
  close();
  if (open()) return 1;
  
  // Calculate the length of each coordinate frame in bytes
  frame_lines = (P->natom * 3) / 10;
  if (((P->natom * 3) % 10) > 0)
    frame_lines++;
  // If DOS, CR present for each newline
  if (File->isDos) frame_lines*=2;
  frameSize = ((P->natom * 3 * 8) + frame_lines);
  // REMD header size is 42 bytes if present
  frameSize += hasREMD;
  if (debug>0) {
    mprintf("Each frame has %i chars plus %i newlines",
            P->natom*3*8,frame_lines);
    if (hasREMD>0) mprintf(" and REMD header");
    mprintf(". Total %i.\n",frameSize);
  }
  // Allocate memory to buffer 1 frame
  frameBuffer=(char*) calloc(frameSize, sizeof(char));
  // Read the first frame of coordinates
  if ( File->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    rprintf("ERROR in read of Coords frame 1\n");
    return 1;
  }
  // DEBUG - Print first line of coords
  if (debug>0) {
    memcpy(buffer,frameBuffer,81);
    buffer[80]='\0';
    rprintf("First 80 bytes: %s\n",buffer);
  }
  /* Check for box coordinates. If present, update the frame size and
   * reallocate the frame buffer.
   */
  if ( File->IO->Gets(buffer,BUFFER_SIZE)==0 ) {
    if (debug>0) rprintf("DEBUG: Line after first frame: (%s)\n",buffer);
    lineSize=strlen(buffer);

    if (strncmp(buffer,"REMD",4)==0 || strncmp(buffer,"HREMD",5)==0) {
      // REMD header - no box coords
      BoxType=0;
    } else if (lineSize<80) {
      /* Line is shorter than 80 chars, indicates box coords.
       * Length of the line HAS to be a multiple of 8, and probably could be
       * enforced to only 3 or 6.
       * Subtract 1 char for newline.
       */
      if (debug>0) mprintf("    Box line is %i chars.\n",lineSize);
      if ( ((lineSize-1)%24)!=0 ) {
        mprintf("Error in box coord line. Expect only 3 or 6 box coords.\n");
        return 1;
      }
      numBoxCoords=(lineSize-1) / 8;
      if (debug>0) mprintf("    Detected %i box coords.\n",numBoxCoords);
      // Determine box type based on angles. Angles are usually not printed 
      // for orthogonal and truncated octahedral boxes, but check here just
      // to be safe. If no angles present use parmtop Box Type.
      if (numBoxCoords>3) {
        sscanf(buffer, "%8lf%8lf%8lf%8lf%8lf%8lf",box,box+1,box+2,box+3,box+4,box+5);
        CheckBoxType(box);
      } else
        BoxType = P->BoxType; 
      // Reallocate frame buffer accordingly
      frameSize+=lineSize;
      frameBuffer=(char*) realloc(frameBuffer,frameSize * sizeof(char));
    }
  }

  /* Calculate number of frames. If not possible and this is not a
   * compressed file the trajectory is probably corrupted. Frames will
   * be read until EOF.
   * NOTE: It is necessary to use the stat command to get the file size
   * instead of fseek in case the file has been popen'd.
   * NOTE: Need the uncompressed file size!
   */
  if (File->uncompressed_size>0)
    file_size = File->uncompressed_size;
  else
    file_size=File->file_size;
  if (debug>0)
    rprintf("Title offset=%i FrameSize=%i UncompressedFileSize=%lli\n",
            titleSize,frameSize,file_size);
  frame_size = (long long int) titleSize;
  file_size = file_size - frame_size; // Subtract title size from file total size.
  frame_size = (long long int) frameSize;
  Frames = (int) (file_size / frame_size);

  if ( (file_size % frame_size) == 0 ) {
    seekable = 1;
    stop = Frames;
  } else {
    stop = -1;
    Frames = -1;
    seekable = 0;
    mprintf("  Error: Could not predict # frames in %s. This usually indicates \n",File->filename);
    mprintf("         a corrupted trajectory. Frames will be read until EOF.\n");
  }

  if (debug>0)
    rprintf("Atoms: %i FrameSize: %i TitleSize: %i NumBox: %i Seekable: %i Frames: %i\n\n", 
            P->natom, frameSize, titleSize, numBoxCoords, seekable, Frames);
  // Close the file
  close();

  return 0;
}
#undef ROUTINE

/*
 * AmberTraj::WriteArgs()
 * Check for remdtraj flag.
 */
int AmberTraj::WriteArgs(ArgList *A) {
  if (A->hasKey("remdtraj")) hasREMD=REMD_HEADER_SIZE;
  return 0;
}

/*
 * AmberTraj::SetupWrite()
 * Set up trajectory for write. Calculate the length of each cooordinate
 * frame. Set the title and titleSize. Calculate the full output file
 * size, necessary only for seeking when MPI writing. Allocate memory for
 * the frame buffer. 
 */
int AmberTraj::SetupWrite() {
  int frame_lines;
  long int outfilesize;

  // Calculate the length of each coordinate frame in bytes
  frame_lines = (P->natom * 3) / 10;
  if (((P->natom * 3) % 10) > 0)
    frame_lines++;
  frameSize = ((P->natom * 3 * 8) + frame_lines);
  // REMD header size is 42 bytes
  frameSize += hasREMD;

  // If box coords are present, allocate extra space for them
  // NOTE: Currently only writing box lengths for all box types. This means
  //       writing triclinic box type is currently not supported.
  if (BoxType!=0) {
    numBoxCoords=3; // Only write out box lengths for trajectories
    frameSize+=((numBoxCoords*8)+1);
  }

  // Set up title and size.  
  if (title==NULL) {
    title=(char*) malloc(81*sizeof(char));
    strcpy(title,"Cpptraj Generated trajectory");
    titleSize=81;
  } else {
    // Check title length - Add one byte for newline char
    titleSize=strlen(title) + 1;
    if (titleSize>81) {
      rprintf("Error: AmberTraj::open: Given title is too long!\n[%s]\n",title);
      return 1;
    }
  }

  // Calculate total output file size.
  // NOTE: Need to correct for instances where # frames not known
  outfilesize=(long int) (titleSize+(P->parmFrames * frameSize));

  File->OpenFile();
  // Set file size - relevant for MPI only
  File->IO->SetSize( outfilesize ); 
  // Write title - master only
  //File->IO->Rank_printf(0,"%-80s\n",title);
  File->CloseFile();

  if (debug>0) {
    rprintf("(%s): Each frame has %i chars plus %i newlines. Total %i. Filesize %li\n", 
            File->filename,P->natom*3*8,frame_lines,frameSize,outfilesize);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer=(char*) calloc(frameSize+1, sizeof(char));

  return 0;
}

/*
 * Info()
 */
void AmberTraj::Info() {
  if (hasREMD) 
    mprintf("is an AMBER REMD trajectory");
  else
    mprintf("is an AMBER trajectory");
}
