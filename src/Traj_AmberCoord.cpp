// Traj_AmberCoord
#include <cstdio>
#include "Traj_AmberCoord.h"
#include "CpptrajStdio.h"

/// Size of REMD header
const size_t Traj_AmberCoord::REMD_HEADER_SIZE = 42;
/// Max Size of 1 line of amber coords, 80 + newline [+CR] + NULL
const size_t Traj_AmberCoord::BUF_SIZE = 128;

// CONSTRUCTOR
Traj_AmberCoord::Traj_AmberCoord() :
  debug_(0),
  natom3_(0),
  titleSize_(0),
  hasREMD_(0),
  numBoxCoords_(0),
  outfmt_("%8.3lf"),
  highPrecision_(false),
  mdvel_(0),
  seekable_(false)
{
  boxAngle_[0] = 0.0; 
  boxAngle_[1] = 0.0; 
  boxAngle_[2] = 0.0; 
}

// DESTRUCTOR
Traj_AmberCoord::~Traj_AmberCoord() {
  if (mdvel_!=0) delete mdvel_;
}

// IsRemdHeader()
static inline bool IsRemdHeader(const char* buffer) {
  if ( (buffer[0]=='R' && buffer[1]=='E' && buffer[2]=='M' && buffer[3]=='D') ||
       (buffer[0]=='H' && buffer[1]=='R' && buffer[2]=='E' && buffer[3]=='M'))
    return true;
  return false;
}

// Traj_AmberCoord::ID_TrajFormat()
bool Traj_AmberCoord::ID_TrajFormat(CpptrajFile& fileIn) {
  char buffer2[BUF_SIZE];
  // File must already be set up for read
  if (fileIn.OpenFile()) return false;
  fileIn.Gets(buffer2, BUF_SIZE); // Title
  fileIn.Gets(buffer2, BUF_SIZE); // REMD header/coords
  fileIn.CloseFile();
  // Check if second line contains REMD/HREMD, Amber Traj with REMD header
  if ( IsRemdHeader( buffer2 ) ) {
    if (debug_>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
    hasREMD_ = REMD_HEADER_SIZE;
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

// Traj_AmberCoord::closeTraj()
void Traj_AmberCoord::closeTraj() {
  file_.CloseFile();
}

// Traj_AmberCoord::openTraj()
/** Open Amber Trajectory, read past title and set titleSize
  * Always read title so file pointer is always correctly positioned
  * to read a frame.
  */
int Traj_AmberCoord::openTraj() {
  char titleIn[BUF_SIZE];
  switch (file_.Access()) {
    case CpptrajFile::READ: // Read in title, set size in bytes 
      if ( file_.OpenFile() ) return 1;
      if ( file_.Gets(titleIn,BUF_SIZE) ) {
        rprintf("Warning: EOF encountered during reading of title from (%s)\n", 
                file_.BaseFileStr());
        return 1;
      }
      title_.assign( titleIn );
      titleSize_ = title_.size();
      break;
    case CpptrajFile::WRITE: // Write title
      if (file_.OpenFile()) return 1;
      file_.Rank_printf(0,"%-*s\n",titleSize_-1,title_.c_str());
      break;
    case CpptrajFile::APPEND:  // Just open
      if (file_.OpenFile()) return 1;
      break;
  }
  return 0;
}

// Traj_AmberCoord::readFrame()
/** Read coordinates from Amber trajectory into Frame. */ 
// NOTE: Precalculate the header, coord, and box offsets.
// NOTE: There are currently NO checks for null for X, box, and T!
int Traj_AmberCoord::readFrame(int set, double *X, double *V, double *box, double *T) {
#ifdef TRAJDEBUG
  mprinterr("Reading frame %10i from %s\n",set,BaseFileStr());
#endif
  // If trajectory is not broken, seek to given frame.
  if (seekable_)
    file_.SeekToFrame( set ); 

  // Read frame into the char buffer
  if (file_.ReadFrame() == -1) return 1;

  // Get REMD Temperature if present
  if (hasREMD_ > 0) 
    file_.GetDoubleAtPosition(*T, 33, 41); 
  // Get Coordinates; offset is hasREMD (size in bytes of REMD header)
  file_.BufferBeginOffset();
  file_.BufferToDouble(X, natom3_);
  if (numBoxCoords_ != 0) { 
    file_.BufferToDouble(box, numBoxCoords_);
    // Set box angles to parmtop default if not read in
    if (numBoxCoords_==3) {
      box[3] = boxAngle_[0];
      box[4] = boxAngle_[1];
      box[5] = boxAngle_[2];
    }
  }
  return 0;
}

// Traj_AmberCoord::writeFrame()
/** Write coordinates from Frame to frameBuffer. frameBuffer must be large
  * enough to accomodate all coords in F (handled by SetupWrite).
  */
// NOTE: The output frame size is calcd here - should it just be precalcd?
int Traj_AmberCoord::writeFrame(int set, double *X, double *V, double *box, double T) {
  if (hasREMD_>0) 
    file_.Printf("REMD  %8i %8i %8i %8.3f\n", 0, set+1, set+1, T);

  file_.BufferBegin();
  file_.DoubleToBuffer(X, natom3_, outfmt_);
  if (numBoxCoords_ != 0) 
    file_.DoubleToBuffer(box, numBoxCoords_, outfmt_);

  if (file_.WriteFrame()) return 1;

  return 0;
}

// Traj_AmberCoord::processReadArgs()
int Traj_AmberCoord::processReadArgs(ArgList& argIn) {
  ArgList::ConstArg mdvelname = argIn.getKeyString("mdvel");
  if (mdvelname!=NULL) {
    // Set up corresponding velocity file
    if (mdvel_!=0) delete mdvel_;
    mdvel_ = new Traj_AmberCoord();
    CpptrajFile velfile;
    if (velfile.SetupRead( mdvelname, debug_ )) {
      mprinterr("Error: mdvel %s: Could not set up file for reading.\n", mdvelname);
      return 1;
    }
    if ( !mdvel_->ID_TrajFormat( velfile ) ) {
      mprinterr("Error: mdvel %s is not Amber ASCII velocity file.\n", mdvelname);
      return 1;
    }
    mprintf("\tmdvel: %s\n", mdvelname);
  }
    
  return 0;
}

// Traj_AmberCoord::setupTrajin()
/** Setup opens the given file for read access, sets information.
  * \return the number of frames present in trajectory file.
  * \return -1 if an error occurs.
  * \return -2 if the number of frames could not be determined.
  */
int Traj_AmberCoord::setupTrajin(std::string const& fname, Topology* trajParm,
                    TrajInfo& tinfo)
{
  // Set up file for reading 
  if (file_.SetupRead( fname, debug_ )) return -1;
  // Attempt to open the file. open() sets the title and titleSize
  if (openTraj()) return -1;
  // Allocate mem to read in frame (plus REMD header if present) 
  natom3_ = trajParm->Natom() * 3;
  file_.SetupFrameBuffer( natom3_, 8, 10, hasREMD_ );
  if (debug_ > 0) {
    mprintf("Each frame is %u bytes", file_.FrameSize());
    if (hasREMD_ > 0) mprintf(" (including REMD header)");
    mprintf(".\n");
  }
  // Read the first frame of coordinates
  if ( file_.ReadFrame() == -1 ) {
    mprinterr("Error in read of Coords frame 1 of trajectory %s.\n", file_.BaseFileStr());
    return -1;
  }
  // Check for box coordinates. If present, update the frame size and
  // reallocate the frame buffer.
  tinfo.BoxInfo.SetNoBox();
  std::string nextLine = file_.GetLineUnbuffered();
  if ( !nextLine.empty() ) {
    if (debug_>0) rprintf("DEBUG: Line after first frame: (%s)\n", nextLine.c_str());
    if ( IsRemdHeader(nextLine.c_str()) ) {
      // REMD header - no box coords
      numBoxCoords_ = 0;
    } else {
      double box[6];
      numBoxCoords_ = sscanf(nextLine.c_str(), "%8lf%8lf%8lf%8lf%8lf%8lf",
                             box, box+1, box+2, box+3, box+4, box+5);
      if (numBoxCoords_ == -1) {
        mprinterr("Error: Could not read Box coord line of trajectory %s.",file_.BaseFileStr());
        return -1;
      } else if (numBoxCoords_ == 3) {
        // Box lengths only, orthogonal or truncated oct. Set lengths.
        tinfo.BoxInfo.SetLengths( box );
      } else if (numBoxCoords_ == 6) {
        // General triclinic. Set lengths and angles.
        tinfo.BoxInfo.SetBox( box );
      } else {
        mprinterr("Error: Expect only 3 or 6 box coords, got %i\n", numBoxCoords_);
        return -1;
      }
    }
    // Reallocate frame buffer accordingly
    file_.ResizeBuffer( numBoxCoords_ );
  }
  // Calculate Frames and determine seekable. If not possible and this is not a
  // compressed file the trajectory is probably corrupted. Frames will
  // be read until EOF (Frames = -2).
  int Frames = 0;
  if (debug_>0)
    rprintf("Title offset=%lu FrameSize=%lu UncompressedFileSize=%lu\n",
            titleSize_, file_.FrameSize(), file_.UncompressedSize());
  off_t title_size = (off_t) titleSize_;
  off_t frame_size = (off_t) file_.FrameSize();
  off_t uncompressed_size = file_.UncompressedSize();
  off_t file_size = uncompressed_size - title_size;
  if (file_.Compression() != CpptrajFile::NO_COMPRESSION) {
    // -----==== AMBER TRAJ COMPRESSED ====------
    // If the uncompressed size of compressed file is reported as <= 0,
    // uncompressed size cannot be determined. Read coordinates until
    // EOF.
    if (uncompressed_size <= 0) {
      mprintf("Warning: %s: Uncompressed size of trajectory could not be determined.\n",
              file_.BaseFileStr());
      if (file_.Compression() == CpptrajFile::BZIP2)
        mprintf("         (This is normal for bzipped files)\n");
      mprintf("         Number of frames could not be calculated.\n");
      mprintf("         Frames will be read until EOF.\n");
      Frames = -2;
      seekable_ = false;
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
          off_t tmpfsize = ((file_size * 4) - uncompressed_size) / 4294967296LL;
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
        seekable_ = true;
      } else {
        mprintf("Warning: %s: Number of frames in compressed traj could not be determined.\n",
                file_.BaseFileStr());
        mprintf("         Frames will be read until EOF.\n");
        Frames=-2;
        seekable_=false;
      }
    }
  } else {     
    // ----==== AMBER TRAJ NOT COMPRESSED ====----
    Frames = (int) (file_size / frame_size);
    if ( (file_size % frame_size) == 0 ) {
      seekable_ = true;
    } else {
      mprintf("Warning: %s: Could not accurately predict # frames. This usually \n",
              file_.BaseFileStr());
      mprintf("         indicates a corrupted trajectory. Will attempt to read %i frames.\n",
              Frames);
      seekable_=false;
    }
  }

  if (debug_>0)
    rprintf("Atoms: %i FrameSize: %lu TitleSize: %lu NumBox: %i Seekable: %i Frames: %i\n\n", 
            trajParm->Natom(), frame_size, title_size, numBoxCoords_, (int)seekable_, Frames);
  // Close the file
  closeTraj();
  // Set trajectory info.
  tinfo.NreplicaDim = 0;
  tinfo.BoxInfo.SetAngles( boxAngle_ );
  tinfo.HasV = false;
  tinfo.HasT = (hasREMD_ > 0);
  tinfo.IsSeekable = seekable_;
  tinfo.Title = title_;

  return Frames;
}

// Traj_AmberCoord::processWriteArgs()
int Traj_AmberCoord::processWriteArgs(ArgList& argIn) {
  if (argIn.hasKey("remdtraj")) hasREMD_ = REMD_HEADER_SIZE; 
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
int Traj_AmberCoord::setupTrajout(std::string const& fname, Topology* trajParm,
                     int NframesToWrite, TrajInfo const& tinfo, bool append)
{
  // Set up file for writing
  if (append) {
    // First set up and open READ to read title. If file does not exist
    // this will fail silently and title will not be set.
    if (file_.SetupRead( fname, debug_ ) == 0) {
      if (openTraj()) return 1; // This sets the title and titleSize
      closeTraj();
    }
    file_.SwitchAccess( CpptrajFile::APPEND );
  } else {
    file_.SetupWrite( fname, debug_ );
    title_ = tinfo.Title;
  }   
  // Set Temperature Write 
  if (tinfo.HasT) hasREMD_ = REMD_HEADER_SIZE;
  // Allocate buffer. Will not buffer REMD header.
  natom3_ = trajParm->Natom() * 3;
  file_.SetupFrameBuffer( natom3_, 8, 10, 0 );
  // If box coords are present, allocate extra space for them
  switch (tinfo.BoxInfo.Type()) {
    case Box::NOBOX   : numBoxCoords_ = 0; break;
    case Box::ORTHO   :
    case Box::TRUNCOCT: numBoxCoords_ = 3; break;
    default           : numBoxCoords_ = 6;
  }
  file_.ResizeBuffer( numBoxCoords_ );
 
  if (title_.empty()) {
    // Set up default title
    title_.assign("Cpptraj Generated trajectory");
    title_.resize(80,' ');
    titleSize_ = title_.size();
  } else {
    // Check title length
    titleSize_ = title_.size(); 
    if (titleSize_>80) {
      title_.resize(80);
      mprintf("Warning: Amber traj title for %s too long: truncating.\n[%s]\n",
              file_.BaseFileStr(), title_.c_str());
    }
  }

  if (debug_>0) 
    rprintf("(%s): Each frame has %lu bytes.\n", file_.BaseFileStr(), file_.FrameSize());
  
  return 0;
}

// Traj_AmberCoord::Info()
void Traj_AmberCoord::info() {
  if (hasREMD_) 
    mprintf("is an AMBER REMD trajectory");
  else
    mprintf("is an AMBER trajectory");
  if (highPrecision_) mprintf(" (high precision)");
}
