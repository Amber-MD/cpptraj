// AmberCoord
#include <cstdlib> // atof
#include <cstdio>  // sprintf
#include <cstring> // strlen
#include "Traj_AmberCoord.h"
#include "CpptrajStdio.h"
// Size of REMD header
const size_t AmberCoord::REMD_HEADER_SIZE = 42;
// Max Size of 1 line of amber coords, 80 + newline [+CR] + NULL
const size_t AmberCoord::BUF_SIZE = 128;

// CONSTRUCTOR
AmberCoord::AmberCoord() {
  //fprintf(stderr,"AmberCoord Constructor.\n");
  titleSize_ = 0;
  frameSize_ = 0;
  hasREMD_ = 0;
  numBoxCoords_ = 0;
  frameBuffer_ = NULL;
  outfmt_ = "%8.3lf";
  highPrecision_ = false;
}

// IsRemdHeader()
static inline bool IsRemdHeader(char* buffer) {
  if ( (buffer[0]=='R' && buffer[1]=='E' && buffer[2]=='M' && buffer[3]=='D') ||
       (buffer[0]=='H' && buffer[1]=='R' && buffer[2]=='E' && buffer[3]=='M'))
    return true;
  return false;
}

bool AmberCoord::ID_TrajFormat() {
  char buffer2[BUF_SIZE];
  // File must already be set up for read
  if (OpenFile()) return false;
  IO->Gets(buffer2, BUF_SIZE); // Title
  IO->Gets(buffer2, BUF_SIZE); // REMD header/coords
  CloseFile();
  // Check if second line contains REMD/HREMD, Amber Traj with REMD header
  if ( IsRemdHeader( buffer2 ) ) {
    if (debug_>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
    hasREMD_ = REMD_HEADER_SIZE;
    hasTemperature_ = true;
    return true;
  }
  // Check if we can read at least 3 coords of width 8, Amber trajectory
  float TrajCoord[3];
  if ( sscanf(buffer2, "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 )
  {
    if (debug_>0) mprintf("  AMBER TRAJECTORY file\n");
    return true;
  }
  return false;
}

// AmberCoord::closeTraj()
void AmberCoord::closeTraj() {
  CloseFile();
}

// AmberCoord::openTraj()
/** Open Amber Trajectory, read past title and set titleSize
  * Always read title so file pointer is always correctly positioned
  * to read a frame.
  */
int AmberCoord::openTraj() {
  char titleIn[BUF_SIZE];
  switch (access_) {
    case READ:
      if (OpenFile()) return 1;
      // Read in title, set size in bytes 
      if ( IO->Gets(titleIn,BUF_SIZE) ) {
        rprintf( "Warning: EOF encountered during reading of title from (%s)\n", BaseName());
        return 1;
      }
      title_.assign( titleIn );
      titleSize_ = title_.size();
      break;
    case WRITE:
      if (OpenFile()) return 1;
      // Write title
      // NOTE: THIS WILL SCREW UP IF OPEN CALLED MULTIPLE TIMES FOR WRITE!
      Rank_printf(0,"%-*s\n",titleSize_-1,title_.c_str());
      break;
    case APPEND:
      // First, open the traj with read access to set the title
      if (title_.empty()) {
        access_ = READ;
        if (openTraj()) return 1;
        closeTraj();
      }
      // Reset the access and open append
      access_ = APPEND;
      if (OpenFile()) return 1;
      break;
  }
  return 0;
}

// AmberCoord::readFrame()
/** Read coordinates from Amber trajectory into Frame. */ 
// NOTE: Precalculate the header, coord, and box offsets.
// NOTE: There are currently NO checks for null for X, box, and T!
int AmberCoord::readFrame(int set, double *X, double *V, double *box, double *T) {
  off_t offset;

#ifdef TRAJDEBUG
  mprinterr("Reading frame %10i from %s\n",set,BaseName());
#endif
 
  // If trajectory is not broken, seek to given frame.
  if (seekable_) {
    offset = (off_t) set;
    offset *= (off_t) frameSize_;
    offset += (off_t) titleSize_;
    IO->Seek(offset);
  }

  // Read frame into the char buffer
  if ( IO->Read(frameBuffer_,frameSize_)==-1 ) return 1;

  // Get REMD Temperature if present
  if (hasREMD_>0) {
    char savechar = frameBuffer_[41];
    frameBuffer_[41] = '\0';
    *T = atof(frameBuffer_ + 33);
    frameBuffer_[41] = savechar;
    //if (debug>0) fprintf(stderr,"DEBUG: Replica T is %lf (%s)\n",F->T,Temp);
  }
  // Get Coordinates - hasREMD is size in bytes of REMD header
  BufferBegin(hasREMD_);
  BufferToDouble(X, natom3_, 8);
  if (hasBox_) { 
    BufferToDouble(box, numBoxCoords_ ,8);
    // Set box angles to parmtop default if not read in
    if (numBoxCoords_==3) {
      box[3] = boxAngle_[0];
      box[4] = boxAngle_[1];
      box[5] = boxAngle_[2];
    }
  }
  return 0;
}

// AmberCoord::SetHighPrecision()
/** Change the output format from 8.3 to 8.6 */
void AmberCoord::SetHighPrecision() {
  outfmt_="%8.6lf";
  highPrecision_=true;
}

// AmberCoord::writeFrame()
/** Write coordinates from Frame to frameBuffer. frameBuffer must be large
  * enough to accomodate all coords in F (handled by SetupWrite).
  */
// NOTE: The output frame size is calcd here - should it just be precalcd?
int AmberCoord::writeFrame(int set, double *X, double *V, double *box, double T) {
  size_t outFrameSize;

  BufferBegin();

  if (hasREMD_>0) {
    sprintf(bufferPosition_,"REMD  %8i %8i %8i %8.3lf\n",0,set+OUTPUTFRAMESHIFT,
            set+OUTPUTFRAMESHIFT,T);
    bufferPosition_ += hasREMD_;
  }

  DoubleToBuffer(X, natom3_, outfmt_, 8, 10);
  if (hasBox_) 
    BoxToBuffer(box, numBoxCoords_, outfmt_, 8);

  outFrameSize = (size_t) (bufferPosition_ - frameBuffer_);
  
  if (IO->Write(frameBuffer_,outFrameSize)) return 1;

  return 0;
}

// AmberCoord::setupTrajin()
/** Setup opens the given file for read access, sets information.
  * \return the number of frames present in trajectory file.
  * \return -1 if an error occurs.
  * \return -2 if the number of frames could not be determined.
  */
int AmberCoord::setupTrajin(Topology *trajParm) {
  char buffer[BUF_SIZE];
  int maxi;
  size_t frame_lines, lineSize;
  off_t file_size, frame_size, tmpfsize, title_size;
  int Frames;
  double box[6]; // For checking box coordinates
  bool sizeFound;

  // Attempt to open the file. open() sets the title and titleSize
  if (openTraj()) return -1;

  // Calculate the length of each coordinate frame in bytes
  natom3_ = trajParm->Natom() * 3;
  frame_lines = (size_t)(natom3_ / 10);
  if ((natom3_ % 10) > 0)
    ++frame_lines;
  // If DOS, CR present for each newline
  if (isDos_) frame_lines *= 2;
  frameSize_ = (((size_t)natom3_ * 8) + frame_lines);
  // REMD header size is 42 bytes if present
  frameSize_ += hasREMD_;
  if (debug_>0) {
    mprintf("Each frame has %i chars plus %lu newlines", natom3_*8,frame_lines);
    if (hasREMD_>0) mprintf(" and REMD header");
    mprintf(". Total %lu.\n",frameSize_);
  }
  // Allocate memory to buffer 1 frame
  frameBuffer_ = new char[ frameSize_ ];
  // Read the first frame of coordinates
  if ( IO->Read(frameBuffer_, frameSize_)==-1 ) {
    mprinterr("Error in read of Coords frame 1 of trajectory %s.\n",BaseName());
    return -1;
  }
  // DEBUG - Print first line of coords
  if (debug_>0) {
    char savechar = frameBuffer_[80];
    frameBuffer_[80]='\0';
    rprintf("First 80 bytes: %s\n",frameBuffer_);
    frameBuffer_[80] = savechar;
  }
  // Check for box coordinates. If present, update the frame size and
  // reallocate the frame buffer.
  if ( IO->Gets(buffer,BUF_SIZE)==0 ) {
    if (debug_>0) rprintf("DEBUG: Line after first frame: (%s)\n",buffer);
    lineSize = strlen(buffer);

    if ( IsRemdHeader(buffer) ) {
      // REMD header - no box coords
      hasBox_ = false;
    } else if (lineSize<80) {
      // Line is shorter than 80 chars, indicates box coords.
      // Length of the line HAS to be a multiple of 8, and probably could be
      // enforced to only 3 or 6. Subtract 1 char for newline.
      if (debug_>0) mprintf("    Box line is %i chars.\n",lineSize);
      if ( ((lineSize-1)%24)!=0 ) {
        mprinterr("Error in box coord line of trajectory %s.\n",BaseName());
        mprinterr("      Expect only 3 or 6 box coords.\n");
        mprinterr("Problem line: %s\n",buffer);
        return -1;
      }
      hasBox_ = true;
      numBoxCoords_ = (int)((lineSize-1) / 8);
      if (debug_>0) mprintf("    Detected %i box coords.\n",numBoxCoords_);
      // If present, get box angles. Angles are usually not printed 
      // for orthogonal and truncated octahedral boxes, but check here just
      // to be safe. These will be used to determine the box type in TrajectoryFile.
      // NOTE: lengths could be an empty read.
      if (numBoxCoords_>3) 
        sscanf(buffer, "%8lf%8lf%8lf%8lf%8lf%8lf",box,box+1,box+2,
               boxAngle_,boxAngle_+1,boxAngle_+2);
      // Reallocate frame buffer accordingly
      frameSize_ += lineSize;
      delete[] frameBuffer_;
      frameBuffer_ = new char[ frameSize_ ];
    }
  }

  // Calculate Frames and determine seekable. If not possible and this is not a
  // compressed file the trajectory is probably corrupted. Frames will
  // be read until EOF (Frames = -2).
  if (debug_>0)
    rprintf("Title offset=%lu FrameSize=%lu Size=%lu UncompressedFileSize=%lu\n",
            titleSize_,frameSize_,file_size_,uncompressed_size_);
  title_size = (off_t) titleSize_;
  frame_size = (off_t) frameSize_;
  // -----==== AMBER TRAJ COMPRESSED ====------
  if (compressType_ != NO_COMPRESSION) {
    // If the uncompressed size of compressed file is reported as <= 0,
    // uncompressed size cannot be determined. Read coordinates until
    // EOF.
    if (uncompressed_size_ <= 0) {
      mprintf("Warning: %s: Uncompressed size of trajectory could not be determined.\n",
              BaseName());
      if (compressType_==BZIP2)
        mprintf("         (This is normal for bzipped files)\n");
      mprintf("         Number of frames could not be calculated.\n");
      mprintf("         Frames will be read until EOF.\n");
      Frames = -2;
      seekable_ = false;
    } else {
      file_size = uncompressed_size_;
      file_size = file_size - title_size;
      // Frame calculation for large gzip files
      if (compressType_ == GZIP) {
        // Since this is gzip compressed, if the file_size % frame size != 0, 
        // it could be that the uncompressed filesize > 4GB. Since 
        //   ISIZE = uncompressed % 2^32, 
        // try ((file_size + (2^32 * i)) % frame_size) and see if any are 0.
        if ( (file_size % frame_size) != 0) {
          // Determine the maximum number of iterations to try based on the
          // fact that Amber trajectories typically compress about 3x with
          // gzip.
          tmpfsize = ((file_size * 4) - uncompressed_size_) / 4294967296LL;
          maxi = (int) tmpfsize;
          ++maxi;
          if (debug_>1)
            mprintf("\tLooking for uncompressed gzip size > 4GB, %i iterations.\n",maxi);
          tmpfsize = 0;
          sizeFound=false;
          for (int i = 0; i < maxi; i++ ) {
            tmpfsize = (4294967296LL * i) + file_size;
            if ( (tmpfsize % frame_size) == 0) {sizeFound=true; break;}
          }
          if (sizeFound) file_size = tmpfsize;
        }
      }
      if ((file_size % frame_size) == 0) {
        Frames = (int) (file_size / frame_size);
        seekable_ = true;
      } else {
        mprintf("Warning: %s: Number of frames in compressed traj could not be determined.\n",
                BaseName());
        mprintf("         Frames will be read until EOF.\n");
        Frames=-2;
        seekable_=false;
      }
    }
  // ----==== AMBER TRAJ NOT COMPRESSED ====----
  } else {     
    file_size = file_size_;
    file_size = file_size - title_size;
    Frames = (int) (file_size / frame_size);
    if ( (file_size % frame_size) == 0 ) {
      seekable_ = true;
    } else {
      mprintf("Warning: %s: Could not accurately predict # frames. This usually \n",
              BaseName());
      mprintf("         indicates a corrupted trajectory. Will attempt to read %i frames.\n",
              Frames);
      seekable_=false;
    }
  }

  if (debug_>0)
    rprintf("Atoms: %i FrameSize: %lu TitleSize: %lu NumBox: %i Seekable: %i Frames: %i\n\n", 
            trajParm->Natom(), frameSize_, titleSize_, numBoxCoords_, (int)seekable_, Frames);
  // Close the file
  closeTraj();

  return Frames;
}

// AmberCoord::processWriteArgs()
int AmberCoord::processWriteArgs(ArgList *argIn) {
  if (argIn->hasKey("remdtraj")) this->SetTemperature();
  if (argIn->hasKey("highprecision")) this->SetHighPrecision();
  return 0;
}

// AmberCoord::setupTrajout()
/** Set up trajectory for write. Calculate the length of each cooordinate
  * frame. Set the title and titleSize. Calculate the full output file
  * size, necessary only for seeking when MPI writing. Allocate memory for
  * the frame buffer. 
  */
int AmberCoord::setupTrajout(Topology *trajParm) {
  int frame_lines;
  //long int outfilesize;

  // Calculate the length of each coordinate frame in bytes
  natom3_ = trajParm->Natom() * 3;
  frame_lines = natom3_ / 10;
  if ((natom3_ % 10) > 0)
    ++frame_lines;
  frameSize_ = (natom3_ * 8) + (size_t)frame_lines;
  // REMD header size is 42 bytes
  if (hasTemperature_) hasREMD_ = REMD_HEADER_SIZE;
  frameSize_ += hasREMD_;

  // If box coords are present, allocate extra space for them
  // NOTE: Currently only writing box lengths for all box types. This means
  //       writing triclinic box type is currently not supported.
  if (hasBox_) {
    numBoxCoords_=3; // Only write out box lengths for trajectories
    frameSize_+=(((size_t)numBoxCoords_*8)+1);
  }

  // Set up title and size.  
  if (title_.empty()) {
    title_.assign("Cpptraj Generated trajectory");
    title_.resize(80,' ');
    titleSize_ = title_.size();
  } else {
    // Check title length
    titleSize_ = title_.size(); 
    if (titleSize_>80) {
      title_.resize(80);
      mprintf("Warning: Amber traj title for %s too long: truncating.\n[%s]\n",
              BaseName(), title_.c_str());
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

  if (debug_>0) {
    //rprintf("(%s): Each frame has %i chars plus %i newlines. Total %i. Filesize %li\n", 
    rprintf("(%s): Each frame has %i chars plus %i newlines. Total %lu.\n", 
            BaseName(),natom3_*8,frame_lines,frameSize_);
  }

  // Allocate memory to buffer 1 frame
  // One extra char for NULL
  frameBuffer_ = new char[ frameSize_ + 1 ];

  return 0;
}

// AmberCoord::Info()
void AmberCoord::info() {
  if (hasREMD_) 
    mprintf("is an AMBER REMD trajectory");
  else
    mprintf("is an AMBER trajectory");
  if (highPrecision_) mprintf(" (high precision)");
}
