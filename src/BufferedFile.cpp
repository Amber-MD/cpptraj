#include <cstdio>  // sprintf
#include <cstdlib> // atof
#include <cstring> // memset
#include "BufferedFile.h"
#include "CpptrajStdio.h"

BufferedFile::BufferedFile() :
  buffer_(0),
  bufferPosition_(0),
  frameSize_(0),
  offset_(0),
  Ncols_(0),
  eltWidth_(0)
{ }

BufferedFile::~BufferedFile() {
  if (buffer_!=0) delete[] buffer_;
}

/** Prepare the buffer to receive organized chunks of text, i.e. 
  * organized in some regular fashion (e.g. an Amber Traj, which
  * is 10 cols of 8.3 precision floating point numbers etc).
  * \param Nelts Total expected number of elements to read.
  * \param eltWidth Width in chars of each element.
  * \param eltsPerLine Number of elements per line (columns).
  * \param offsetIn Offset to be used in seeking.
  * \return Size of set-up frame.
  */
size_t BufferedFile::SetupFrameBuffer(int Nelts, int eltWidthIn, int eltsPerLine, int offsetIn) 
{
  Ncols_ = eltsPerLine;
  eltWidth_ = (size_t)eltWidthIn;
  offset_ = (size_t) offsetIn;
  int frame_lines = Nelts / eltsPerLine;
  if ((Nelts % eltsPerLine) > 0)
    ++frame_lines;
  bool readingFile = (Access() == CpptrajFile::READ);
  // If Reading and DOS, CR present for each newline
  if (readingFile && IsDos()) frame_lines *= 2;
  // Calculate total frame size
  frameSize_ = (((size_t)Nelts * eltWidth_) + frame_lines) + offset_;
  // If writing, add +1 for NULL
  if (!readingFile)
    ++frameSize_;
  if (buffer_!=0) delete[] buffer_;
  if (frameSize_ < 1) 
    buffer_ = 0;
  else {
    buffer_ = new char[ frameSize_ ];
    memset(buffer_, 0, frameSize_);
  }
  bufferPosition_ = buffer_;
  return frameSize_;
}

/** Increase size of buffer by delta elements, keeping contents
  * intact. Intended for adding one line (e.g. for box coords).
  */
size_t BufferedFile::ResizeBuffer(int delta) {
  if (delta == 0) return frameSize_;
  if (delta < 0) {
    mprinterr("Internal Error: ResizeBuffer: Negative value given.\n");
    return frameSize_;
  }
  // Add 1 for newline
  size_t newsize = frameSize_ + ((size_t)delta * eltWidth_) + 1;
  // If reading and isDos add 1 for carriage return
  if (Access() == CpptrajFile::READ && IsDos())
    ++newsize;
  char* newbuffer = new char[ newsize ];
  memcpy(newbuffer, buffer_, frameSize_);
  memset(newbuffer+frameSize_, 0, newsize - frameSize_);
  delete[] buffer_;
  buffer_ = newbuffer;
  bufferPosition_ = buffer_;
  frameSize_ = newsize;
  return frameSize_;
}

int BufferedFile::SeekToFrame(size_t set) {
  return Seek( (off_t)((set * frameSize_) + offset_) );
}

int BufferedFile::ReadFrame() {
  return Read( buffer_, frameSize_ );
}

int BufferedFile::WriteFrame() {
  return Write( buffer_, (size_t)(bufferPosition_ - buffer_) );
}

void BufferedFile::GetDoubleAtPosition(double& val, size_t start, size_t end) {
  char savechar = buffer_[end];
  buffer_[end] = '\0';
  val = atof(buffer_ + start);
  buffer_[end] = savechar;
}

void BufferedFile::BufferBeginOffset() {
  bufferPosition_ = buffer_ + offset_;
}

void BufferedFile::BufferBegin() {
  bufferPosition_ = buffer_;
}

void BufferedFile::BufferBeginAt(size_t pos) {
  bufferPosition_ = buffer_ + pos;
}

/** Convert text in buffer containing numerical elements with format 
  * X0Y0Z0X1Y1Z1...XNYNZN to the given double array. The width of each 
  * element should be what SetupFrameBuffer was called with, and the 
  * number of elements to read should not be greater than Nelts.
  * Newlines are skipped. Output array should be as big as Nout. 
  * Update bufferPosition after read.
  */
void BufferedFile::BufferToDouble(double* Xout, int Nout) {
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

/** Convert given double array to ordered text in buffer. The number of
  * elements in the given array should be what SetupFrameBuffer was
  * called with. Update bufferPosition after write. 
  */
void BufferedFile::DoubleToBuffer(const double* Xin, int Nin, const char* format)
{
  int col = 0;
  for (int element = 0; element < Nin; ++element) {
    sprintf(bufferPosition_, format, Xin[element]);
    bufferPosition_ += eltWidth_;
    ++col;
    if ( col == Ncols_ ) {
      sprintf(bufferPosition_,"\n");
      ++bufferPosition_;
      col = 0;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( col != 0 ) {
    sprintf(bufferPosition_,"\n");
    ++bufferPosition_;
  }
}

