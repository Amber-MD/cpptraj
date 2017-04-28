// Traj_AmberCoord
#include <cstdio> // sscanf
#include "Traj_AmberCoord.h"
#include "CpptrajStdio.h"

/// Size of REMD header
const size_t Traj_AmberCoord::REMD_HEADER_SIZE = 42;
const size_t Traj_AmberCoord::RXSGLD_HEADER_SIZE = 44;

// CONSTRUCTOR
Traj_AmberCoord::Traj_AmberCoord() :
  natom3_(0),
  headerSize_(0),
  numBoxCoords_(0),
  outfmt_("%8.3lf"),
  highPrecision_(false),
  outputTemp_(false)
{
  boxAngle_[0] = 0.0; 
  boxAngle_[1] = 0.0; 
  boxAngle_[2] = 0.0; 
}

// IsRemdHeader()
static inline bool IsRemdHeader(const char* buffer) {
  if ( (buffer[0]=='R' && buffer[1]=='E' && buffer[2]=='M' && buffer[3]=='D') ||
       (buffer[0]=='H' && buffer[1]=='R' && buffer[2]=='E' && buffer[3]=='M'))
    return true;
  return false;
}

static inline bool IsRxsgldHeader(const char* buffer) {
  if ( buffer[0]=='R' && buffer[1]=='X' && buffer[2]=='S' && buffer[3]=='G' )
    return true;
  return false;
}

// Traj_AmberCoord::ID_TrajFormat()
bool Traj_AmberCoord::ID_TrajFormat(CpptrajFile& fileIn) {
  // File must already be set up for read
  if (fileIn.OpenFile()) return false;
  if (fileIn.NextLine()==0) return false; // Title
  std::string buffer2 = fileIn.GetLine(); // REMD header/coords
  fileIn.CloseFile();
  // Check if second line contains REMD/HREMD, Amber Traj with REMD header
  if ( IsRemdHeader( buffer2.c_str() ) ) {
    if (debug_>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
    headerSize_ = REMD_HEADER_SIZE + (size_t)fileIn.IsDos();
    tStart_ = 33; // 42 - 8 - 1
    tEnd_   = 41; // 42 - 1
    return true;
  }
  // TODO: Read these in as indices instead of temperatures
  if ( IsRxsgldHeader( buffer2.c_str() ) ) {
    mprintf("  AMBER TRAJECTORY with RXSGLD header.\n");
    headerSize_ = RXSGLD_HEADER_SIZE + (size_t)fileIn.IsDos();
    tStart_ = 35; // 44 - 8 - 1
    tEnd_   = 43; // 44 - 1
    return true;
  }
  // Check if we can read at least 3 coords of width 8, Amber trajectory
  float TrajCoord[3];
  if ( sscanf(buffer2.c_str(), "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 )
  {
    if (debug_>0) mprintf("  AMBER TRAJECTORY file\n");
    return true;
  }
  return false;
}

// Traj_AmberCoord::closeTraj()
void Traj_AmberCoord::closeTraj() {
  file_.CloseFile();
}

// Traj_AmberCoord::openTrajin()
/** Open Amber Trajectory for reading. Read past title; this is done so
  * file pointer is always correctly positioned to read first frame.
  */
int Traj_AmberCoord::openTrajin() {
  if (file_.OpenFile()) return 1;
  if (file_.NextLine() == 0) return 1; // Skip past title
  return 0;
}

// Traj_AmberCoord::readFrame()
/** Read coordinates from Amber trajectory into Frame. */ 
// NOTE: Precalculate the header, coord, and box offsets.
// NOTE: There are currently NO checks for null for X, box, and T!
int Traj_AmberCoord::readFrame(int set, Frame& frameIn) {
#ifdef TRAJDEBUG
  mprinterr("Reading frame %10i from %s\n",set,Filename().base());
#endif
  // Seek to given frame.
  file_.SeekToFrame( set ); 

  // Read frame into the char buffer
  if (file_.ReadFrame()) return 1;

  // Get REMD Temperature if present
  if (headerSize_ != 0) 
    file_.GetDoubleAtPosition(*(frameIn.tAddress()), tStart_, tEnd_); 
  // Get Coordinates; offset is hasREMD (size in bytes of REMD header)
  file_.BufferBeginAt(headerSize_);
  file_.BufferToDouble(frameIn.xAddress(), natom3_);
  if (numBoxCoords_ != 0) { 
    file_.BufferToDouble(frameIn.bAddress(), numBoxCoords_);
    // Set box angles to parmtop default if not read in
    if (numBoxCoords_==3)
      frameIn.SetBoxAngles( boxAngle_ );
  }
  return 0;
}

int Traj_AmberCoord::readVelocity(int set, Frame& frameIn) {
  file_.SeekToFrame( set );
  // Read frame into the char buffer
  if (file_.ReadFrame()) return 1;
  file_.BufferBeginAt(headerSize_);
  file_.BufferToDouble(frameIn.vAddress(), natom3_);
  return 0;
}

// Traj_AmberCoord::writeFrame()
/** Write coordinates from Frame to frameBuffer. frameBuffer must be large
  * enough to accomodate all coords in F (handled by SetupWrite).
  */
// NOTE: The output frame size is calcd here - should it just be precalcd?
int Traj_AmberCoord::writeFrame(int set, Frame const& frameOut) {
  if (headerSize_ != 0) 
    file_.Printf("REMD  %8i %8i %8i %8.3f\n", 0, set+1, set+1, frameOut.Temperature());

  file_.BufferBegin();
  file_.DoubleToBuffer(frameOut.xAddress(), natom3_, outfmt_);
  if (numBoxCoords_ != 0) 
    file_.DoubleToBuffer(frameOut.bAddress(), numBoxCoords_, outfmt_);

  if (file_.WriteFrame()) return 1;

  return 0;
}

// Traj_AmberCoord::setupTrajin()
int Traj_AmberCoord::setupTrajin(FileName const& fname, Topology* trajParm)
{
  // Set up file for reading 
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  // Attempt to open the file. Read and set the title and titleSize
  if ( file_.OpenFile() ) return TRAJIN_ERR;
  std::string title = file_.GetLine();
  // Allocate mem to read in frame (plus REMD header if present). REMD
  // header is checked for when file is IDd. Title size is used in seeking. 
  natom3_ = trajParm->Natom() * 3;
  file_.SetupFrameBuffer( natom3_, 8, 10, headerSize_, title.size() );
  if (debug_ > 0) {
    mprintf("Each frame is %u bytes", file_.FrameSize());
    if (headerSize_ != 0) mprintf(" (including REMD header)");
    mprintf(".\n");
  }
  // Read the first frame of coordinates
  if ( file_.ReadFrame() ) {
    mprinterr("Error: in read of Coords frame 1 of trajectory %s.\n", file_.Filename().base());
    return TRAJIN_ERR;
  }
  // Check for box coordinates. If present, update the frame size and
  // reallocate the frame buffer. If less than 3 atoms there is no way
  // to tell if a line is a box line or coordinates, so skip.
  Box boxInfo;
  if ( trajParm->Natom() < 3 ) {
    mprintf("Warning: Less than 3 atoms, skipping box check.\n");
    numBoxCoords_ = 0;
  } else {
    std::string nextLine = file_.GetLine();
    if ( !nextLine.empty() ) {
      if (debug_>0) rprintf("DEBUG: Line after first frame: (%s)\n", nextLine.c_str());
      if ( IsRemdHeader(nextLine.c_str()) || IsRxsgldHeader(nextLine.c_str()) ) {
        // REMD header - no box coords
        numBoxCoords_ = 0;
      } else {
        double box[8];
        numBoxCoords_ = sscanf(nextLine.c_str(), "%8lf%8lf%8lf%8lf%8lf%8lf%8lf%8lf",
                               box, box+1, box+2, box+3, box+4, box+5, box+6, box+7);
        if (numBoxCoords_ == -1) {
          mprinterr("Error: Could not read Box coord line of trajectory %s.\n",
                   file_.Filename().base());
          return TRAJIN_ERR;
        } else if (numBoxCoords_ == 8) {
          // Full line of coords was read, no box coords.
          numBoxCoords_ = 0;
        } else if (numBoxCoords_ == 3) {
          // Box lengths only, ortho. or truncated oct. Use default parm angles.
          if (trajParm->ParmBox().Type() == Box::NOBOX)
            mprintf("Warning: Trajectory only contains box lengths and topology has no box info.\n"
                    "Warning: To set box angles for topology use the 'parmbox' command.\n");
          box[3] = boxAngle_[0] = trajParm->ParmBox().Alpha();
          box[4] = boxAngle_[1] = trajParm->ParmBox().Beta();
          box[5] = boxAngle_[2] = trajParm->ParmBox().Gamma();
          boxInfo.SetBox( box );
        } else if (numBoxCoords_ == 6) {
          // General triclinic. Set lengths and angles.
          boxInfo.SetBox( box );
        } else {
          mprinterr("Error: In %s, expect only 3 or 6 box coords, got %i\n"
                    "Error:   Box line=[%s]\n",
                    file_.Filename().base(), numBoxCoords_, nextLine.c_str());
          return TRAJIN_ERR;
        }
      }
    }
    // Reallocate frame buffer accordingly
    file_.ResizeBuffer( numBoxCoords_ );
  }
  // Calculate Frames and determine seekable. If not possible and this is not a
  // compressed file the trajectory is probably corrupted. Frames will
  // be read until EOF.
  int Frames = 0;
  if (debug_>0)
    rprintf("Title offset=%lu FrameSize=%lu UncompressedFileSize=%lu\n",
            title.size(), file_.FrameSize(), file_.UncompressedSize());
  off_t title_size = (off_t) title.size();
  off_t frame_size = (off_t) file_.FrameSize();
  off_t uncompressed_size = file_.UncompressedSize();
  off_t file_size = uncompressed_size - title_size;
  bool seekable = false;
  if (file_.Compression() != CpptrajFile::NO_COMPRESSION) {
    // -----==== AMBER TRAJ COMPRESSED ====------
    // If the uncompressed size of compressed file is reported as <= 0,
    // uncompressed size cannot be determined. Read coordinates until
    // EOF.
    if (uncompressed_size <= 0) {
      mprintf("Warning: %s: Uncompressed size of trajectory could not be determined.\n",
              file_.Filename().base());
      if (file_.Compression() == CpptrajFile::BZIP2)
        mprintf("         (This is normal for bzipped files)\n");
      mprintf("         Number of frames could not be calculated.\n");
      mprintf("         Frames will be read until EOF.\n");
      Frames = TRAJIN_UNK;
      seekable = false;
    } else {
      // Frame calculation for large gzip files
      if (file_.Compression() == CpptrajFile::GZIP) {
        // Since this is gzip compressed, if the file_size % frame size != 0, 
        // it could be that the uncompressed filesize > 4GB. Since 
        //   ISIZE = uncompressed % 2^32, 
        // try ((file_size + (2^32 * i)) % frame_size) and see if any are 0.
        if ( (file_size % frame_size) != 0) {
          // Determine the maximum number of iterations to try based on the
          // fact that Amber trajectories typically compress about 3x with
          // gzip.
          off_t tmpfsize = ((file_.FileSize() * 4) - uncompressed_size) / 4294967296LL;
          int maxi = (int) tmpfsize;
          ++maxi;
          if (debug_>1)
            mprintf("\tLooking for uncompressed gzip size > 4GB, %i iterations.\n",maxi);
          tmpfsize = 0;
          bool sizeFound = false;
          for (int i = 0; i < maxi; i++ ) {
            tmpfsize = (4294967296LL * i) + file_size;
            if ( (tmpfsize % frame_size) == 0) {sizeFound=true; break;}
          }
          if (sizeFound) file_size = tmpfsize;
        }
      }
      if ((file_size % frame_size) == 0) {
        Frames = (int) (file_size / frame_size);
        seekable = true;
      } else {
        mprintf("Warning: %s: Number of frames in compressed traj could not be determined.\n"
                "Warning:  Frames will be read until EOF.\n", file_.Filename().base());
        Frames = TRAJIN_UNK;
        seekable = false;
      }
    }
  } else {     
    // ----==== AMBER TRAJ NOT COMPRESSED ====----
    Frames = (int) (file_size / frame_size);
    if ( (file_size % frame_size) == 0 ) {
      seekable = true;
    } else {
      mprintf("Warning: %s: Could not accurately predict # frames. This usually\n"
              "Warning:  indicates a corrupted trajectory or trajectory/topology\n"
              "Warning:  mismatch. Will attempt to read %i frames.\n",
              file_.Filename().base(), Frames);
      seekable = false;
    }
  }

  if (debug_>0)
    rprintf("Atoms: %i FrameSize: %lu TitleSize: %lu NumBox: %i Seekable: %i Frames: %i\n\n", 
            trajParm->Natom(), frame_size, title_size, numBoxCoords_, (int)seekable, Frames);
  // Close the file
  file_.CloseFile();
  // Set trajectory info: no velocity, no time.
  SetCoordInfo( CoordinateInfo(boxInfo, false, (headerSize_ != 0), false) );
  SetTitle( title );
  return Frames;
}

void Traj_AmberCoord::WriteHelp() {
  mprintf("\tremdtraj     : Write temperature to trajectory (makes REMD trajectory).\n"
          "\thighprecision: (ADVANCED USE ONLY) Write 8.6 instead of 8.3 format.\n");
}

// Traj_AmberCoord::processWriteArgs()
int Traj_AmberCoord::processWriteArgs(ArgList& argIn) {
  outputTemp_ = argIn.hasKey("remdtraj");
  if (argIn.hasKey("highprecision")) { 
    outfmt_ = "%8.6lf";
    highPrecision_ = true;
  }
  return 0;
}

// Traj_AmberCoord::setupTrajout()
/** Set up trajectory for write. Calculate the length of each cooordinate
  * frame. Set the title and titleSize. Calculate the full output file
  * size, necessary only for seeking when MPI writing. Allocate memory for
  * the frame buffer. 
  */
int Traj_AmberCoord::setupTrajout(FileName const& fname, Topology* trajParm,
                                  CoordinateInfo const& cInfoIn, 
                                  int NframesToWrite, bool append)
{
  // Set Temperature Write
  // FIXME: Check for temperatures in input frames and when appending.
  SetCoordInfo( cInfoIn );
  if (outputTemp_) {
    headerSize_ = REMD_HEADER_SIZE;
    if (!CoordInfo().HasTemp())
      mprintf("Warning: No temperature information in input frames.\n");
  }
  if (!append) {
    // Write the title if not appending
    if (file_.SetupWrite( fname, debug_ )) return 1;
    std::string title = Title();
    if (title.empty()) {
      // Set up default title
      title.assign("Cpptraj Generated trajectory");
      title.resize(80,' ');
      SetTitle( title );
    } else {
      // Check title length
      if (title.size() > 80) {
        mprintf("Warning: Amber traj title for %s too long: truncating.\n[%s]\n",
                file_.Filename().base(), title.c_str());
        title.resize(80);
      }
    }
    if (file_.OpenFile()) return 1;
    file_.Printf("%-s\n", title.c_str());
  } else {
    // Just set up for append
    if (file_.SetupAppend( fname, debug_ )) return 1;
    if (file_.OpenFile()) return 1;
  }
  // Allocate buffer. Will not buffer REMD header or need to seek.
  // NOTE: This is done here since SetupFrameBuffer allocates based
  //       on write mode, which is not known until now.
  natom3_ = trajParm->Natom() * 3;
  file_.SetupFrameBuffer( natom3_, 8, 10 );
  // If box coords are present, allocate extra space for them
  switch (CoordInfo().TrajBox().Type()) {
    case Box::NOBOX   : numBoxCoords_ = 0; break;
    case Box::ORTHO   :
    case Box::TRUNCOCT: numBoxCoords_ = 3; break;
    default           : numBoxCoords_ = 6;
  }
  file_.ResizeBuffer( numBoxCoords_ );
 
  if (debug_>0) 
    rprintf("(%s): Each frame has %lu bytes.\n", file_.Filename().base(), file_.FrameSize());
  
  return 0;
}

// Traj_AmberCoord::Info()
void Traj_AmberCoord::Info() {
  if (CoordInfo().HasTemp()) 
    mprintf("is an AMBER REMD trajectory");
  else
    mprintf("is an AMBER trajectory");
  if (highPrecision_) mprintf(" (high precision)");
}
#ifdef MPI
// =============================================================================
int Traj_AmberCoord::parallelOpenTrajin(Parallel::Comm const& commIn) {
  mprinterr("Error: Parallel read not supported for Amber ASCII coords.\n");
  return 1;
}

/** This assumes file has been previously set up with parallelSetupTrajout
  * and title has been written, so open append.
  */
int Traj_AmberCoord::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return (file_.ParallelOpenFile( CpptrajFile::APPEND, commIn ));
}

/** First master performs all necessary setup, then sends info to all children.
  */
int Traj_AmberCoord::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{
  int err = 0;
  // In parallel MUST know # of frames to write in order to correctly set size
  if (NframesToWrite < 1) {
    mprinterr("Error: # frames to write must be known for Amber Coords output in parallel.\n");
    err = 1;
  } else if (commIn.Master()) {
    // NOTE: This writes the title.
    err = setupTrajout(fname, trajParm, cInfoIn, NframesToWrite, append);
    // NOTE: setupTrajout leaves file open. Should this change?
    file_.CloseFile();
  }
  commIn.MasterBcast(&err, 1, MPI_INT);
  if (err != 0) return 1;
  // Synchronize info on non-master threads.
  SyncTrajIO( commIn );
  // TODO For simplicity convert everything to double. Is this just lazy?
  double tmpArray[8];
  if (commIn.Master()) {
    tmpArray[0] = (double)natom3_;
    tmpArray[1] = (double)headerSize_;
    tmpArray[2] = (double)tStart_;
    tmpArray[3] = (double)tEnd_;
    tmpArray[4] = (double)numBoxCoords_;
    tmpArray[5] = boxAngle_[0];
    tmpArray[6] = boxAngle_[1];
    tmpArray[7] = boxAngle_[2];
    commIn.MasterBcast(tmpArray, 8, MPI_DOUBLE);
  } else {
    commIn.MasterBcast(tmpArray, 8, MPI_DOUBLE);
    natom3_       = (int)tmpArray[0];
    headerSize_   = (size_t)tmpArray[1];
    tStart_       = (size_t)tmpArray[2];
    tEnd_         = (size_t)tmpArray[3];
    numBoxCoords_ = (int)tmpArray[4];
    boxAngle_[0]  = tmpArray[5];
    boxAngle_[1]  = tmpArray[6];
    boxAngle_[2]  = tmpArray[7];
    if (append)
      file_.SetupAppend( fname, debug_ );
    else
      file_.SetupWrite( fname, debug_ );
  }
  // For parallel output we will need to seek. Set up the buffer again with correct offsets.
  // Figure out the size of the written title.
  unsigned int titleSize = (unsigned int)Title().size() + 1; // +1 for newline
  titleSize = std::min(81U, titleSize);
  file_.SetupFrameBuffer( natom3_, 8, 10, headerSize_, titleSize );
  file_.ResizeBuffer( numBoxCoords_ );
  if (debug_>0)
    rprintf("'%s'(Parallel): Each frame has %lu bytes.\n", file_.Filename().base(),
            file_.FrameSize());
  // TODO set file size
  return 0;
}

int Traj_AmberCoord::parallelReadFrame(int set, Frame& frameIn) { return 1; }

int Traj_AmberCoord::parallelWriteFrame(int set, Frame const& frameOut) {
  // Seek to given frame.
  file_.SeekToFrame( set );
  return ( writeFrame(set, frameOut) );
}

void Traj_AmberCoord::parallelCloseTraj() { closeTraj(); }
#endif
