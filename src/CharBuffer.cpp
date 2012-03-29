// CharBuffer
#include <cstring> // memcpy
#include <cstdio>  // sprintf
#include <cstdlib> // atof
#include <cstdarg> // va_list
#include <cctype>  // isdigit
#include <vector>
#include "CharBuffer.h"
// DEBUG
#include "CpptrajStdio.h"

// CONSTRUCTOR
CharBuffer::CharBuffer() {
  buffer_ = NULL;
  ptr_ = NULL;
  bufferSize_ = 0;
  end_buffer_ = NULL;
}

// DESTRUCTOR
CharBuffer::~CharBuffer() {
  if (buffer_!=NULL) 
    delete[] buffer_;
}

// Copy Constructor
CharBuffer::CharBuffer(const CharBuffer& rhs) {
  bufferSize_ = rhs.bufferSize_;
  buffer_ = new char[ bufferSize_ ];
  memcpy(buffer_, rhs.buffer_, bufferSize_*sizeof(char));
  ptr_ = buffer_ + (rhs.ptr_ - rhs.buffer_);
  end_buffer_ = buffer_ + bufferSize_;
}

// Assignment
CharBuffer &CharBuffer::operator=(const CharBuffer &rhs) {
  if (this == &rhs) return *this;
  if (buffer_!=NULL) delete[] buffer_;
  bufferSize_ = rhs.bufferSize_;
  buffer_ = new char[ bufferSize_ ];
  memcpy(buffer_, rhs.buffer_, bufferSize_*sizeof(char));
  ptr_ = buffer_ + (rhs.ptr_ - rhs.buffer_);
  end_buffer_ = buffer_ + bufferSize_;
  return *this;
}

// CharBuffer::Buffer()
/** Return a pointer to the beginning of buffer.
  */
char *CharBuffer::Buffer() {
  return buffer_;
}

// CharBuffer::c_str()
/** Like Buffer() except ensure a NULL terminated string is returned.
  */
const char *CharBuffer::c_str() {
  *ptr_='\0';
  return buffer_;
}

// CharBuffer::BufferSize()
size_t CharBuffer::BufferSize() { 
  return bufferSize_; 
}

// CharBuffer::BufferPtr() 
char *CharBuffer::BufferPtr() { 
  return ptr_; 
}

// CharBuffer::CurrentSize()
size_t CharBuffer::CurrentSize() {
  //printf("DEBUG:\tCurrentSize=%lu, Allocd for %lu\n",(size_t) (ptr - buffer),bufferSize);
  return (size_t) (ptr_ - buffer_);
}

// CharBuffer::Allocate()
/// Allocate space for a character buffer of size <sizeIn>.
void CharBuffer::Allocate(size_t sizeIn) {
  if (buffer_!=NULL) delete[] buffer_;
  bufferSize_ = sizeIn;
  buffer_ = new char[ bufferSize_ ];
  ptr_ = buffer_;
  end_buffer_ = buffer_ + bufferSize_;
}

// CharBuffer::IncreaseSize()
/** Increase the size of the char buffer by delta, ensuring that the existing
  * data is kept. Place the ptr at the beginning of the new region.
  * NOTE: If used after Sprintf this can result in uninitialized bytes.
  */
void CharBuffer::IncreaseSize(size_t delta) {
  size_t newsize = bufferSize_ + delta;
  char *newbuffer = new char[ newsize ];
  if (buffer_!=NULL) {
    memcpy(newbuffer, buffer_, bufferSize_*sizeof(char));
    delete[] buffer_;
  }
  buffer_ = newbuffer;
  ptr_ = buffer_ + bufferSize_; 
  bufferSize_ = newsize;
  end_buffer_ = buffer_ + bufferSize_;
}

// CharBuffer::Reallocate()
/** Increase size of the char buffer by delta, ensuring existing data is
  * kept. Do not alter position of ptr.
  */
void CharBuffer::Reallocate(size_t delta) {
  size_t newsize = bufferSize_ + delta;
  size_t old_ptr = ptr_ - buffer_; 
  char *newbuffer = new char[ newsize ];
  if (buffer_!=NULL) {
    memcpy(newbuffer, buffer_, bufferSize_*sizeof(char));
    delete[] buffer_;
  }
  buffer_ = newbuffer;
  ptr_ = buffer_ + old_ptr;
  bufferSize_ = newsize;
  end_buffer_ = buffer_ + bufferSize_;
}

// CharBuffer::DumpBuffer()
/// For debugging only, write buffer contents to stdout.
void CharBuffer::DumpBuffer() {
  mprintf("Allocated size=%lu\n",bufferSize_);
  mprintf("Current Size=%lu\n",CurrentSize());
  for (unsigned int i = 0; i < CurrentSize(); i++)
    mprintf("\t%i [%c]\n",i,buffer_[i]);
}

// CharBuffer::Rewind()
void CharBuffer::Rewind() {
  ptr_ = buffer_;
}

// CharBuffer::WriteDouble()
/// Write a double to the buffer with given format.
void CharBuffer::WriteDouble(const char *format, double dval) {
  //printf("DEBUG:\tWriteDouble[%s, %lf]\n",format, dval);
  int n_char_written = sprintf(ptr_, format, dval);
  ptr_ += n_char_written;
}

// CharBuffer::WriteInteger()
/// Write an integer to the buffer with given format.
void CharBuffer::WriteInteger(const char *format, int ival) {
  //printf("DEBUG:\tWriteInteger[%s, %i]\n",format, ival);
  int n_char_written = sprintf(ptr_, format, ival);
  ptr_ += n_char_written;
}

// CharBuffer::WriteString()
/// Write a string to the buffer with given format.
void CharBuffer::WriteString(const char *format, const char *sval) {
  //printf("DEBUG:\tWriteString[%s, %s]\n",format, sval);
  int n_char_written = sprintf(ptr_, format, sval);
  ptr_ += n_char_written;
}

// CharBuffer::WriteXY()
void CharBuffer::WriteXY(const char *format, double X, double Y) {
  int n_char_written = sprintf(ptr_, format, X, Y);
  ptr_ += n_char_written;
}

// CharBuffer::WriteDoubleXYZ()
/// Write the given double array of size 3 to the buffer with given format.
void CharBuffer::WriteDoubleXYZ(const char *format, double *XYZ) {
  int n_char_written = sprintf(ptr_, format, XYZ[0], XYZ[1], XYZ[2]);
  ptr_ += n_char_written;
}

// CharBuffer::NewLine()
/// Add a newline character to buffer.
void CharBuffer::NewLine() {
  ptr_[0]='\n';
  ++ptr_;
}

// CharBuffer::Space()
/// Add space to buffer
void CharBuffer::Space() {
  ptr_[0]=' ';
  ++ptr_;
}

// CharBuffer::Sprintf()
/** Write a string to buffer with printf like syntax. This routine will
  * dynamically reallocate the buffer as needed.
  */
void CharBuffer::Sprintf(const char *fmt, ... ) {
  va_list vl;
  size_t n_chars_left = (size_t)(end_buffer_ - ptr_);
  // Check that buffer wont be blown
  va_start(vl, fmt);
  size_t n_char_written = (size_t)vsnprintf(ptr_, n_chars_left, fmt, vl);
  if (n_char_written<0 || n_char_written >= n_chars_left) {
    Reallocate( n_char_written+1 );
    // Try again
    n_chars_left = (size_t)(end_buffer_ - ptr_);
    va_end(vl);
    va_start(vl, fmt);
    n_char_written = (size_t)vsnprintf(ptr_, n_chars_left, fmt, vl);
  }
  va_end(vl);
  ptr_ += n_char_written;
}  

// CharBuffer::Read()
int CharBuffer::Read(void *str, size_t numbytes) {
  size_t bytes_to_read;
  // Figure out actual size to copy
  size_t bytes_remaining = end_buffer_ - ptr_;
  if (bytes_remaining < 1)
    return -1;
  else if (numbytes > bytes_remaining)
    bytes_to_read = bytes_remaining;
  else
    bytes_to_read = numbytes;

  // Read from ptr into str
  memcpy(str, ptr_, bytes_to_read);
  // Advance ptr
  ptr_ += bytes_to_read;
  return (int)bytes_to_read;
}
 
// CharBuffer::Gets()
/// Get specified number of bytes from the buffer, up to end or NULL.
// [D][o][g][\n][\0]
//  0  1  2  3   4
// String length: 4, array size: 5
// [D][o][g][\0]
//  0  1  2  3
// String length: 3, array size: 4 
int CharBuffer::Gets(char *str, int num) {
  int currentNum = 0;
  while (ptr_ < end_buffer_ && currentNum < num) {
    str[currentNum++] = *ptr_;
    if (*ptr_=='\n' || *ptr_=='\0') {++ptr_; break;}
    ++ptr_;
  }
  if (currentNum==0) return 1;
  if (currentNum==num) --currentNum;
  str[currentNum]='\0';
  return 0;
}

// =============================================================================
// BufferToDouble()
/** Store character buffer containing XYZ coords with format 
  * X0Y0Z0X1Y1Z1...XNYNZN to the given double array. Each coord has given width,
  * newlines are skipped. buffer should be as big as N x width chars. 
  * Return position in buffer after read. If '*' encountered this indicates
  * overflow in trajectory, return NULL.
  */
/*char *BufferToDouble(char *buffer, double *X, int N, int width) {
  char *ptrbegin = buffer;
  char *ptrend = buffer;
  for (int atom = 0; atom < N; atom++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    // NOTE: Search for '*'??
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    X[atom] = atof(ptrbegin);
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return ptrbegin;
}*/

// DoubleToBuffer()
/** Given an array of double, format, and character width corresponding
  * to format, write coords in array to buffer. 
  * Return the position in the buffer after write. 
  */
/*char *DoubleToBuffer(char *buffer, double *X, int N, const char *format, 
                                 int width, int numCols) 
{
  int coord;
  char *ptr = buffer;
  for (coord=0; coord<N; coord++) {
    sprintf(ptr,format,X[coord]);
    ptr+=width;
    if ( ((coord+1)%numCols)==0 ) {
      sprintf(ptr,"\n");
      ptr++;
    }
  }
  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(ptr,"\n");
    ptr++;
  }

  // Calculate frame size
  //coord = (int) (ptr - buffer);
  //return coord;
  return ptr;
}*/

// BoxToBuffer()
/** Given an array of double[6], format, and character width corresponding
  * to format, write box coords to buffer.
  * Return the position in the buffer after write.
  */
/*char *BoxToBuffer(char *buffer, double *box, int numBox, 
                              const char *format, int width) {
//  int coord;
  char *ptr;

  ptr=buffer;
  // Box
  sprintf(ptr,format,box[0]); ptr+=width;
  sprintf(ptr,format,box[1]); ptr+=width;
  sprintf(ptr,format,box[2]); ptr+=width;
  if (numBox>3) {
    sprintf(ptr,format,box[3]); ptr+=width;
    sprintf(ptr,format,box[4]); ptr+=width;
    sprintf(ptr,format,box[5]); ptr+=width;
  }
  sprintf(ptr,"\n");
  ptr++;

  // Calculate frame size
  //coord = (int) (ptr - buffer);
  //return coord;
  return ptr;
}
*/
