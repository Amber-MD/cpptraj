// AmberCoord
#include "Traj_AmberCoord.h"
#include <cstdlib>
#include <cstring>
#include <cstdio> // sprintf
#include "CpptrajStdio.h"
#include "CharBuffer.h"
// Size of REMD header
#define REMD_HEADER_SIZE 42

// CONSTRUCTOR
AmberCoord::AmberCoord() {
  //fprintf(stderr,"AmberCoord Constructor.\n");
  titleSize=0;
  frameSize=0;
  hasREMD=0;
  numBoxCoords=0;
}

// DESTRUCTOR
AmberCoord::~AmberCoord() {
  //fprintf(stderr,"AmberCoord Destructor.\n");
}    

/* AmberCoord::closeTraj()
 */ 
void AmberCoord::closeTraj() {
  tfile->CloseFile();
}

/* AmberCoord::openTraj()
 * Open Amber Trajectory, read past title and set titleSize
 * Always read title so file pointer is always correctly positioned
 * to read a frame.
 */
int AmberCoord::openTraj() {

  switch (tfile->access) {
    case READ:
      if (tfile->OpenFile()) return 1;
      // Read in title, set size in bytes 
      if (title==NULL)  
        title=(char*) malloc(BUFFER_SIZE * sizeof(char));
      if ( tfile->IO->Gets(title,BUFFER_SIZE) ) {
        rprintf( "Warning: EOF encountered during reading of title from (%s)\n", tfile->filename);
        return 1;
      }
      titleSize=strlen(title);
      break;

    case WRITE:
      if (tfile->OpenFile()) return 1;
      // Write title
      // NOTE: THIS WILL SCREW UP IF OPEN CALLED MULTIPLE TIMES FOR WRITE!
      tfile->IO->Rank_printf(0,"%-*s\n",titleSize-1,title);
      break;

    case APPEND:
      // First, open the traj with read access to set the title
      if (title==NULL) {
        tfile->access=READ;
        if (openTraj()) return 1;
        closeTraj();
      }
      // Reset the access and open append
      tfile->access=APPEND;
      tfile->OpenFile();
      break;
  }

  return 0;
}

/* AmberCoord::readFrame()
 * Read coordinates from Amber trajectory into Frame. 
 * NOTE: Precalculate the header, coord, and box offsets.
 * NOTE: There are currently NO checks for null for X, box, and T!
 */
int AmberCoord::readFrame(int set, double *X, double *box, double *T) {
  char Temp[9];
  char *bufferPosition;
  off_t offset;
 
  // If trajectory is not broken, seek to given frame.
  if (seekable) {
    offset = (off_t) set;
    offset *= (off_t) frameSize;
    offset += (off_t) titleSize;
    tfile->IO->Seek(offset);
  }

  // Read frame into the char buffer
  if ( tfile->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) return 1;

  // Get REMD Temperature if present
  if (hasREMD>0) {
    strncpy(Temp,frameBuffer+33,8);
    Temp[8]='\0';
    *T = atof(Temp);
    //if (debug>0) fprintf(stderr,"DEBUG: Replica T is %lf (%s)\n",F->T,Temp);
  }
  // Get Coordinates - hasREMD is size in bytes of REMD header
  if ( (bufferPosition = BufferToDouble(frameBuffer+hasREMD,X,natom3,8))==NULL ) {
    mprinterr("Error: AmberCoord::readFrame: * detected in coordinates of %s\n",tfile->filename);
    return 1;  
  } 
  if (hasBox) { 
    if ( (bufferPosition = BufferToDouble(bufferPosition,box,numBoxCoords,8))==NULL ) {
      mprinterr("Error: AmberCoord::readFrame: * detected in box coordinates of %s\n",
              tfile->filename);
      return 1;
    }
    // Set box angles to parmtop default if not read in
    if (numBoxCoords==3) {
      box[3] = boxAngle[0];
      box[4] = boxAngle[1];
      box[5] = boxAngle[2];
    }
  }
  return 0;
}

/* AmberCoord::writeFrame()
 * Write coordinates from Frame to frameBuffer. frameBuffer must be large
 * enough to accomodate all coords in F (handled by SetupWrite).
 * NOTE: The output frame size is calcd here - should it just be precalcd?
 */
int AmberCoord::writeFrame(int set, double *X, double *box, double T) {
  int outFrameSize;
  char *bufferPosition;
  //off_t offset;

  bufferPosition = frameBuffer;

  if (hasREMD) {
    sprintf(bufferPosition,"REMD  %8i %8i %8i %8.3lf\n",0,set+OUTPUTFRAMESHIFT,
            set+OUTPUTFRAMESHIFT,T);
    bufferPosition += hasREMD;
  }

  bufferPosition = DoubleToBuffer(bufferPosition,X,natom3,"%8.3lf",8,10);
  if (hasBox) 
    bufferPosition = BoxToBuffer(bufferPosition,box,numBoxCoords,"%8.3lf",8);

  outFrameSize = (int) (bufferPosition - frameBuffer);
  
  //if (seekable) 
  // NOTE: Seek only needs to happen when traj file changes
  //offset = (off_t) currentFrame;
  //offset *= (off_t) outFrameSize;
  //offset += (off_t) titleSize;
  //tfile->IO->Seek( offset);

  if (tfile->IO->Write(frameBuffer,sizeof(char),outFrameSize)) return 1;

  //currentFrame++;
 
  return 0;
}

/* AmberCoord::setupRead()
 * Setup opens the given file for read access, sets information.
 * Return the number of frames present in trajectory file.
 * Return -1 if an error occurs, -2 if the number of frames could 
 * not be determined.
 */
int AmberCoord::setupRead(int natom) {
  char buffer[BUFFER_SIZE];
  int frame_lines;
  int lineSize;
  long long int file_size, frame_size;
  int Frames;
  double box[6]; // For checking box coordinates

  // Attempt to open the file. open() sets the title and titleSize
  if (openTraj()) return -1;

  // Read 1 line to check for REMD header
  if ( tfile->IO->Gets(buffer, BUFFER_SIZE) ) {
    mprinterr("Error: Could not reading second line of trajectory %s.\n",tfile->filename);
    return -1;
  }
  if (strncmp(buffer,"REMD",4)==0 || strncmp(buffer,"HREMD",5)==0) {
    hasREMD=REMD_HEADER_SIZE;
    hasTemperature=true;
  }
  // Reopen the file
  closeTraj();
  if (openTraj()) return -1;
  
  // Calculate the length of each coordinate frame in bytes
  natom3 = natom * 3;
  frame_lines = natom3 / 10;
  if ((natom3 % 10) > 0)
    frame_lines++;
  // If DOS, CR present for each newline
  if (tfile->isDos) frame_lines*=2;
  frameSize = ((natom3 * 8) + frame_lines);
  // REMD header size is 42 bytes if present
  frameSize += hasREMD;
  if (debug>0) {
    mprintf("Each frame has %i chars plus %i newlines", natom3*8,frame_lines);
    if (hasREMD>0) mprintf(" and REMD header");
    mprintf(". Total %i.\n",frameSize);
  }
  // Allocate memory to buffer 1 frame
  frameBuffer=(char*) calloc(frameSize, sizeof(char));
  // Read the first frame of coordinates
  if ( tfile->IO->Read(frameBuffer,sizeof(char),frameSize)==-1 ) {
    mprinterr("Error in read of Coords frame 1 of trajectory %s.\n",tfile->filename);
    return -1;
  }
  // DEBUG - Print first line of coords
  if (debug>0) {
    memcpy(buffer,frameBuffer,81);
    buffer[80]='\0';
    rprintf("First 80 bytes: %s\n",buffer);
  }
  // Check for box coordinates. If present, update the frame size and
  // reallocate the frame buffer.
  if ( tfile->IO->Gets(buffer,BUFFER_SIZE)==0 ) {
    if (debug>0) rprintf("DEBUG: Line after first frame: (%s)\n",buffer);
    lineSize=strlen(buffer);

    if (strncmp(buffer,"REMD",4)==0 || strncmp(buffer,"HREMD",5)==0) {
      // REMD header - no box coords
      hasBox = false;
    } else if (lineSize<80) {
      // Line is shorter than 80 chars, indicates box coords.
      // Length of the line HAS to be a multiple of 8, and probably could be
      // enforced to only 3 or 6. Subtract 1 char for newline.
      if (debug>0) mprintf("    Box line is %i chars.\n",lineSize);
      if ( ((lineSize-1)%24)!=0 ) {
        mprinterr("Error in box coord line of trajectory %s.\n",tfile->filename);
        mprinterr("      Expect only 3 or 6 box coords.\n");
        mprinterr("Problem line: %s\n",buffer);
        return -1;
      }
      hasBox = true;
      numBoxCoords=(lineSize-1) / 8;
      if (debug>0) mprintf("    Detected %i box coords.\n",numBoxCoords);
      // If present, get box angles. Angles are usually not printed 
      // for orthogonal and truncated octahedral boxes, but check here just
      // to be safe. These will be used to determine the box type in TrajectoryFile.
      // NOTE: lengths could be an empty read.
      if (numBoxCoords>3) 
        sscanf(buffer, "%8lf%8lf%8lf%8lf%8lf%8lf",box,box+1,box+2,boxAngle,boxAngle+1,boxAngle+2);
      // Reallocate frame buffer accordingly
      frameSize+=lineSize;
      frameBuffer=(char*) realloc(frameBuffer,frameSize * sizeof(char));
    }
  }

  // Calculate number of frames. If not possible and this is not a
  // compressed file the trajectory is probably corrupted. Frames will
  // be read until EOF.
  // NOTE: It is necessary to use the stat command to get the file size
  // instead of fseek in case the file has been popen'd.
  // NOTE: Need the uncompressed file size!
  if (tfile->uncompressed_size>0)
    file_size = tfile->uncompressed_size;
  else
    file_size=tfile->file_size;
  if (debug>0)
    rprintf("Title offset=%i FrameSize=%i UncompressedFileSize=%lli\n",
            titleSize,frameSize,file_size);
  frame_size = (long long int) titleSize;
  file_size = file_size - frame_size; // Subtract title size from file total size.
  frame_size = (long long int) frameSize;
  Frames = (int) (file_size / frame_size);

  if ( (file_size % frame_size) == 0 ) {
    seekable = true;
  } else {
    Frames = -2;
    seekable = false;
    //mprintf("  Error: Could not predict # frames in %s. This usually indicates \n",tfile->filename);
    //mprintf("         a corrupted trajectory. Frames will be read until EOF.\n");
  }

  if (debug>0)
    rprintf("Atoms: %i FrameSize: %i TitleSize: %i NumBox: %i Seekable: %i Frames: %i\n\n", 
            natom, frameSize, titleSize, numBoxCoords, (int)seekable, Frames);
  // Close the file
  closeTraj();

  return Frames;
}

/* AmberCoord::SetremdTraj()
 * Set hasREMD to REMD_HEADER_SIZE. Triggers write of temperature.
 */
void AmberCoord::SetRemdTraj() {
  hasREMD=REMD_HEADER_SIZE;
}

/* AmberCoord::setupWrite()
 * Set up trajectory for write. Calculate the length of each cooordinate
 * frame. Set the title and titleSize. Calculate the full output file
 * size, necessary only for seeking when MPI writing. Allocate memory for
 * the frame buffer. 
 */
int AmberCoord::setupWrite(int natom) {
  int frame_lines;
  //long int outfilesize;

  // Calculate the length of each coordinate frame in bytes
  natom3 = natom * 3;
  frame_lines = natom3 / 10;
  if ((natom3 % 10) > 0)
    frame_lines++;
  frameSize = (natom3 * 8) + frame_lines;
  // REMD header size is 42 bytes
  frameSize += hasREMD;

  // If box coords are present, allocate extra space for them
  // NOTE: Currently only writing box lengths for all box types. This means
  //       writing triclinic box type is currently not supported.
  if (hasBox) {
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
      rprintf("Error: AmberCoord::open: Given title is too long!\n[%s]\n",title);
      return 1;
    }
  }

  // Calculate total output file size.
  // NOTE: Need to correct for instances where # frames not known
  //outfilesize=(long int) (titleSize+(P->parmFrames * frameSize));

  //tfile->OpenFile();
  // Set file size - relevant for MPI only
  //tfile->IO->SetSize( outfilesize ); 
  // Write title - master only
  //tfile->IO->Rank_printf(0,"%-80s\n",title);
  //tfile->CloseFile();

  if (debug>0) {
    //rprintf("(%s): Each frame has %i chars plus %i newlines. Total %i. Filesize %li\n", 
    rprintf("(%s): Each frame has %i chars plus %i newlines. Total %i.\n", 
            tfile->filename,natom3*8,frame_lines,frameSize);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer=(char*) calloc(frameSize+1, sizeof(char));

  return 0;
}

/* AmberCoord::Info()
 */
void AmberCoord::info() {
  if (hasREMD) 
    mprintf("is an AMBER REMD trajectory");
  else
    mprintf("is an AMBER trajectory");
}
