#include "FrameBuffer.h"
#include <cstring> // memcpy
#include <cstdio> // sprintf
#include <cstdlib> // atof

// CONSTRUCTOR
FrameBuffer::FrameBuffer() :
  frameBuffer_(0),
  bufferPosition_(0),
  frameSize_(0)
{}

// DESTRUCTOR
FrameBuffer::~FrameBuffer() {
  if (frameBuffer_!=0) delete[] frameBuffer_;
}

// Copy Constructor
FrameBuffer::FrameBuffer(const FrameBuffer& rhs) :
  frameBuffer_(0),
  bufferPosition_(0),
  frameSize_(0)
{
  if (rhs.frameBuffer_ != 0) {
    frameSize_ = rhs.frameSize_;
    frameBuffer_ = new char[ frameSize_ ];
    memcpy(frameBuffer_, rhs.frameBuffer_, frameSize_ * sizeof(char));
    bufferPosition_ = frameBuffer_ + (rhs.bufferPosition_ - rhs.frameBuffer_);
  }
}

// Assignment
FrameBuffer& FrameBuffer::operator=(const FrameBuffer& rhs) {
  if (&rhs == this) return *this;
  if (frameBuffer_!=0) delete[] frameBuffer_;
  frameBuffer_ = 0;
  bufferPosition_ = 0;
  frameSize_ = 0;
  if (rhs.frameBuffer_ != 0) {
    frameSize_ = rhs.frameSize_;
    frameBuffer_ = new char[ frameSize_ ];
    memcpy(frameBuffer_, rhs.frameBuffer_, frameSize_ * sizeof(char));
    bufferPosition_ = frameBuffer_ + (rhs.bufferPosition_ - rhs.frameBuffer_);
  }
  return *this;
}

// FrameBuffer::BufferBegin()
void FrameBuffer::BufferBegin(size_t offset) {
  bufferPosition_ = frameBuffer_ + offset;
}

void FrameBuffer::BufferBegin() {
  bufferPosition_ = frameBuffer_;
}

// FrameBuffer::BufferToDouble()
/** Store frameBuffer containing XYZ coords with format 
  * X0Y0Z0X1Y1Z1...XNYNZN to the given double array. Each coord has given width,
  * newlines are skipped. buffer should be as big as N x width chars. Update
  * bufferPosition after read.
  */
// TODO: If '*' encountered indicate error?.
void FrameBuffer::BufferToDouble(double *X, int N, int width) {
  for (int atom = 0; atom < N; atom++) {
    // Advance past newlines / CR (dos)
    while (*bufferPosition_=='\n' || *bufferPosition_=='\r')
      ++bufferPosition_;
    // NOTE: Search for '*'??
    char *ptrend = bufferPosition_ + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    X[atom] = atof(bufferPosition_);
    *ptrend = lastchar;
    bufferPosition_ = ptrend;
  }
}

// FrameBuffer::DoubleToBuffer()
/** Given an array of double, format, and character width corresponding
  * to format, write coords in array to frameBuffer. 
  */
void FrameBuffer::DoubleToBuffer(double *X, int N, const char *format,
                                 int width, int numCols)
{
  int coord = 0;
  for (; coord<N; coord++) {
    sprintf(bufferPosition_,format,X[coord]);
    bufferPosition_ += width;
    if ( ((coord+1)%numCols)==0 ) {
      sprintf(bufferPosition_,"\n");
      ++bufferPosition_;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(bufferPosition_,"\n");
    ++bufferPosition_;
  }

  // Calculate frame size
  //coord = (int) (bufferPosition_ - frameBuffer);
  //return coord;
}

// FrameBuffer::BoxToBuffer()
/** Given an array of double[6], format, and character width corresponding
  * to format, write box coords to frameBuffer.
  */
void FrameBuffer::BoxToBuffer(double *box, int numBox,
                              const char *format, int width) {
  // Box
  sprintf(bufferPosition_,format,box[0]);
  bufferPosition_+=width;
  sprintf(bufferPosition_,format,box[1]);
  bufferPosition_+=width;
  sprintf(bufferPosition_,format,box[2]);
  bufferPosition_+=width;
  if (numBox>3) {
    sprintf(bufferPosition_,format,box[3]);
    bufferPosition_+=width;
    sprintf(bufferPosition_,format,box[4]);
    bufferPosition_+=width;
    sprintf(bufferPosition_,format,box[5]);
    bufferPosition_+=width;
  }
  sprintf(bufferPosition_,"\n");
  ++bufferPosition_;

  // Calculate frame size
  //coord = (int) (bufferPosition_ - frameBuffer);
  //return coord;
}
  
