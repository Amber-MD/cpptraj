#include <cstdio>  // snprintf
#include <cstdlib> // atof
#include "BufferedFrame.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
BufferedFrame::BufferedFrame() :
  buffer_(0),
  bufferPosition_(0),
  frameSize_(0),
  headerSize_(0),
  offset_(0),
  memSize_(0),
  maxSize_(0),
  Ncols_(0),
  col_(0),
  eltWidth_(0),
  saveChar_(0),
  errorCount_(0),
  overflowCount_(0)
{}

/** DESTRUCTOR */
BufferedFrame::~BufferedFrame() {
  if (buffer_!=0) delete[] buffer_;
}

// BufferedFrame::SetupFrameBuffer()
size_t BufferedFrame::SetupFrameBuffer(int Nelts, int eltWidthIn, int eltsPerLine)
{
  return SetupFrameBuffer(Nelts, eltWidthIn, eltsPerLine, 0, 0);
}

// BufferedFrame::SetupFrameBuffer()
size_t BufferedFrame::SetupFrameBuffer(int Nelts, TextFormat const& fmtIn, int eltsPerLine)
{
  writeFmt_ = fmtIn;
  return SetupFrameBuffer(Nelts, writeFmt_.Width(), eltsPerLine, 0, 0);
}
  
/** Prepare the buffer to receive organized chunks of text, i.e. 
  * organized in some regular fashion (e.g. an Amber Traj, which
  * is 10 cols of 8.3 precision floating point numbers etc).
  * \param Nelts Total expected number of elements to read.
  * \param eltWidth Width in chars of each element.
  * \param eltsPerLine Number of elements per line (columns).
  * \param headerSizeIn Any additional bytes at the beginning of frame.
  * \param offsetIn Number of bytes at beginning of file to skip (not part of any frame).
  * \return Size of set-up frame.
  */
size_t BufferedFrame::SetupFrameBuffer(int Nelts, int eltWidthIn, int eltsPerLine, 
                                      size_t headerSizeIn, int offsetIn)
{
  //if (Access() != CpptrajFile::READ &&
  //    buffer_ != 0 && bufferPosition_ != 0 && buffer_ != bufferPosition_)
  //  mprinterr("DEBUG: Buffer was not flushed.\n"); // DEBUG warning
  Ncols_ = eltsPerLine;
  eltWidth_ = (size_t)eltWidthIn;
  offset_ = (size_t) offsetIn;
  headerSize_ = headerSizeIn;
  frameSize_ = CalcFrameSize( Nelts );
  memSize_ = frameSize_ + 1; // +1 for null, TODO not necessary for read?
  //mprintf("DEBUG: Buffer required size %zu, max size %zu.\n", memSize_, maxSize_);
  if (memSize_ > maxSize_) {
    //mprintf("DEBUG: Reallocating.\n");
    // Need to reallocate
    if (buffer_ != 0) delete[] buffer_;
    buffer_ = new char[ memSize_ ];
    maxSize_ = memSize_;
  }
  // Initialize buffer
  std::fill(buffer_, buffer_ + memSize_, 0);
  bufferPosition_ = buffer_;
  col_ = 0;
  saveChar_ = 0;
  //mprintf("DEBUG: %s %i cols, eltWidth= %zu, offset= %zu, frameSize= %zu additional= %zu\n",
  //        Filename().base(), Ncols_, eltWidth_, offset_, frameSize_, headerSize_);
  return frameSize_;
}

/** Based on the current values of Ncols and eltWidth, calculate size
  * in bytes necessary to buffer Nelts.
  */
size_t BufferedFrame::CalcFrameSize( int Nelts ) const {
  int frame_lines = Nelts / Ncols_;
  if ((Nelts % Ncols_) > 0) 
    ++frame_lines;
  bool readingFile = (Access() == CpptrajFile::READ);
  // If Reading and DOS, CR present for each newline
  if (readingFile && IsDos()) frame_lines *= 2;
  // Calculate total frame size
  size_t fsize = (((size_t)Nelts * eltWidth_) + frame_lines) + headerSize_;
  return fsize;
}

/** Increase size of buffer by delta elements, keeping contents intact. */
size_t BufferedFrame::ResizeBuffer(int delta) {
  if (delta == 0) return frameSize_;
  if (delta < 0) {
    mprinterr("Internal Error: ResizeBuffer: Negative value given.\n");
    return frameSize_;
  }
  size_t newsize = frameSize_ + CalcFrameSize( delta );
  char* newbuffer = new char[ newsize + 1]; // +1 for null
  std::copy(buffer_, buffer_+frameSize_, newbuffer);
  std::fill(newbuffer+frameSize_, newbuffer+newsize, 0);
  delete[] buffer_;
  buffer_ = newbuffer;
  bufferPosition_ = buffer_;
  col_ = 0;
  frameSize_ = newsize;
  return frameSize_;
}

// -----------------------------------------------------------------------------
/** Seek to specified frame in file. */
int BufferedFrame::SeekToFrame(size_t set) {
  return Seek( (off_t)((set * frameSize_) + offset_) );
}

/** Set buffer pointer to beginning of buffer. */
void BufferedFrame::BufferBegin() {
  bufferPosition_ = buffer_;
  col_ = 0;
}

/** Set buffer pointer to beginning of buffer, past any header bytes. */
void BufferedFrame::BufferBeginAfterHeader() {
  bufferPosition_ = buffer_ + headerSize_;
  col_ = 0;
}

// -----------------------------------------------------------------------------
/** Attempt to read in the next frameSize_ bytes.
  * \return the actual number of bytes read.
  */
int BufferedFrame::AttemptReadFrame() {
  return Read( buffer_, frameSize_ );
}

/** Read the next frameSize_ bytes.
  * \return true if the number of bytes read was less than frameSize_.
  * \return false if read was successful.
  */
bool BufferedFrame::ReadFrame() {
  int nread = Read(buffer_, frameSize_);
  //mprintf("DEBUG: %i bytes read, frameSize %zu.\n", nread, frameSize_);
  if (nread != (int)frameSize_) {
    if (nread < 0)
      mprinterr("Error: Read %i bytes, expected %zu\n", nread, frameSize_);
    return true;
  }
  return false;
}

/** Convert text in buffer containing numerical elements with format 
  * X0Y0Z0X1Y1Z1...XNYNZN to the given double array. The width of each 
  * element should be what SetupFrameBuffer was called with, and the 
  * number of elements to read should not be greater than Nelts.
  * Newlines are skipped. Output array should be as big as Nout. 
  * Update bufferPosition after read.
  */
void BufferedFrame::BufferToDouble(double* Xout, int Nout) {
  for (int element = 0; element < Nout; ++element) {
    // Advance past newlines / CR (dos)
    while (*bufferPosition_=='\n' || *bufferPosition_=='\r')
      ++bufferPosition_;
    if (*bufferPosition_ == '*') {
      mprinterr("Error: '*' encountered (atom %i", (element / 3) + 1);
      int problem_xyz = element % 3;
      if (problem_xyz == 0)      mprinterr(" X");
      else if (problem_xyz == 1) mprinterr(" Y");
      else                       mprinterr(" Z");
      mprinterr("). This indicates coordinate overflow.\n");
    }
    char *ptrend = bufferPosition_ + eltWidth_;
    char lastchar = *ptrend;
    *ptrend = '\0';
    Xout[element] = atof(bufferPosition_);
    *ptrend = lastchar;
    bufferPosition_ = ptrend;
  }
}

/** \return Pointer to next null-terminated element in buffer.
  */
const char* BufferedFrame::NextElement() {
  if (saveChar_ != 0) *bufferPosition_ = saveChar_;
  const char* position = bufferPosition_;
  bufferPosition_ += eltWidth_;
  //mprinterr("DEBUG: bufferPosition is %zu, frame size is %zu\n", bufferPosition_-buffer_, frameSize_);
  char* end = bufferPosition_;
  while (*bufferPosition_=='\n' || *bufferPosition_=='\r')
    ++bufferPosition_;
  if (bufferPosition_ == end) {
    saveChar_ = *bufferPosition_;
    *bufferPosition_ = '\0';
  } else {
    saveChar_ = 0;
    *end = '\0';
  }
  return position; 
}

// -----------------------------------------------------------------------------
/** Write the current buffer to file. */
int BufferedFrame::WriteFrame() {
  // Report any errors or warnings.
  if (errorCount_ > 0)
    mprinterr("Error: %i errors encountered writing elements to file '%s'\n", errorCount_, Filename().base());
  if (overflowCount_ > 0)
    mprintf("Warning: Character overflow detected writing %i elements to file '%s'\n", overflowCount_, Filename().base());
  errorCount_ = 0;
  overflowCount_ = 0;
  return Write( buffer_, (size_t)(bufferPosition_ - buffer_) );
}

/** Convert given double array to ordered text in buffer. The number of
  * elements in the given array should be what SetupFrameBuffer was
  * called with. Update bufferPosition after write. 
  */
void BufferedFrame::DoubleToBuffer(const double* Xin, int Nin, const char* format)
{
  int col = 0;
  for (int element = 0; element < Nin; ++element) {
    int n_chars = snprintf(bufferPosition_, eltWidth_+1, format, Xin[element]);
    // NOTE: Technically there is an overflow if n_chars is greater than the
    //       element width. However, it is not uncommon for there to be some
    //       basic rounding with floating point formats (particularly in the
    //       mantissa), e.g. '-1013.0374146' (13 chars) becomes
    //       -1013.037415 (12 chars), but snprintf will still return 13.
    //       Therefore, allow for single-digit overflow by checking for
    //       n_chars > eltWidth_+1 instead of n_chars > eltWidth_.
    if (n_chars < 0) {
      errorCount_++;
#     ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
      mprinterr("Error: Writing element '%i' to '%s'\n", element+1, Filename().base());
#     endif
    } else if ((unsigned int)n_chars > eltWidth_+1) {
      overflowCount_++;
#     ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
      mprintf("Warning: Number overflow in '%s', element %i = %f (only writing %zu of %i chars).\n",
              Filename().base(), element+1, Xin[element], eltWidth_, n_chars);
      // DEBUG
      char tmpbuf[128];
      sprintf(tmpbuf, format, Xin[element]);
      mprintf("DEBUG: Full element= '%s'\n", tmpbuf);
#     endif
    }
    bufferPosition_ += eltWidth_;
    ++col;
    if ( col == Ncols_ ) {
      *bufferPosition_ = '\n';
      ++bufferPosition_;
      col = 0;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( col != 0 ) {
    *bufferPosition_ = '\n';
    ++bufferPosition_;
  }
}

/** Advance to the next column of the buffer. Write a newline if 
  * at the end of a line. */
void BufferedFrame::AdvanceCol() {
  bufferPosition_ += eltWidth_;
  ++col_;
  if ( col_ == Ncols_ ) {
    *bufferPosition_ = '\n';
    ++bufferPosition_;
    col_ = 0;
  }
}

/** Write given integer to the buffer and advance to next column. */
void BufferedFrame::IntToBuffer(int ival) {
  int n_chars = snprintf(bufferPosition_, eltWidth_+1, writeFmt_.fmt(), ival);
  if (n_chars < 0) {
    errorCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprinterr("Error: Writing integer %i to '%s'\n", ival, Filename().base());
#   endif
  } else if ((unsigned int)n_chars > eltWidth_) {
    overflowCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprintf("Warning: Number overflow in '%s' integer %i (only writing %zu of %i chars).\n",
            Filename().base(), ival, eltWidth_, n_chars);
#   endif
  }
  AdvanceCol();
}

/** Write the given double to the buffer and advance to next column. */
void BufferedFrame::DblToBuffer(double dval) {
  int n_chars = snprintf(bufferPosition_, eltWidth_+1, writeFmt_.fmt(), dval);
  if (n_chars < 0) {
    errorCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprinterr("Error: Writing double %f to '%s'\n", dval, Filename().base());
#   endif
  } else if ((unsigned int)n_chars > eltWidth_) {
    overflowCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprintf("Warning: Number overflow in '%s' double %g (only writing %zu of %i chars).\n",
            Filename().base(), dval, eltWidth_, n_chars);
#   endif
  }
  AdvanceCol();
}

/** Write the given character string to the buffer and advance to next column. */
void BufferedFrame::CharToBuffer(const char* cval) {
  int n_chars = snprintf(bufferPosition_, eltWidth_+1, writeFmt_.fmt(), cval);
  if (n_chars < 0) {
    errorCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprinterr("Error: Writing string %s to '%s'\n", cval, Filename().base());
#   endif
  } else if ((unsigned int)n_chars > eltWidth_) {
    overflowCount_++;
#   ifdef CPPTRAJ_DEBUG_BUFFEREDFRAME
    mprintf("Warning: Overflow in '%s' string %s (only writing %zu of %i chars).\n",
            Filename().base(), cval, eltWidth_, n_chars);
#   endif
  }
  AdvanceCol();
}

/** Write the buffer to the file. */
void BufferedFrame::FlushBuffer() {
  // If the coord record didnt end on a newline, print one
  if ( col_ != 0 ) {
    *bufferPosition_ = '\n';
    ++bufferPosition_;
  }
  WriteFrame();
  col_ = 0;
  bufferPosition_ = buffer_;
}
