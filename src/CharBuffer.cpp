// CharBuffer
#include <cstring> // memcpy
#include <cstdio>  // sprintf
#include <cstdlib> // atof
#include <cstdarg> // va_list
#include "CharBuffer.h"

// CONSTRUCTOR
CharBuffer::CharBuffer() {
  buffer = NULL;
  ptr = NULL;
  bufferSize = 0;
  end_buffer = NULL;
}

// DESTRUCTOR
CharBuffer::~CharBuffer() {
  if (buffer!=NULL) delete[] buffer;
}

// CharBuffer::Allocate()
/// Allocate space for a character buffer of size <sizeIn>.
void CharBuffer::Allocate(size_t sizeIn) {
  if (buffer!=NULL) delete[] buffer;
  bufferSize = sizeIn;
  buffer = new char[ bufferSize ];
  ptr = buffer;
  end_buffer = buffer + bufferSize;
}

// CharBuffer::IncreaseSize()
/** Increase the size of the char buffer by delta, ensuring that the existing
  * data is kept. Place the ptr at the beginning of the new region.
  */
void CharBuffer::IncreaseSize(size_t delta) {
  size_t newsize = bufferSize + delta;
  char *newbuffer = new char[ newsize ];
  if (buffer!=NULL) {
    memcpy(newbuffer, buffer, bufferSize);
    delete[] buffer;
  }
  buffer = newbuffer;
  ptr = buffer + bufferSize; 
  bufferSize = newsize;
  end_buffer = buffer + bufferSize;
}

// CharBuffer::CurrentSize()
/// Return the size of the data that has been currently written to the buffer.
size_t CharBuffer::CurrentSize() {
  //printf("DEBUG:\tCurrentSize=%lu, Allocd for %lu\n",(size_t) (ptr - buffer),bufferSize);
  return (size_t) (ptr - buffer);
}

// CharBuffer::Rewind()
/// Set pointer to beginning of buffer.
void CharBuffer::Rewind() {
  ptr = buffer;
}

 

// CharBuffer::Sprintf()
/// Write data to the buffer with printf-like syntax.
void CharBuffer::Sprintf(const char *fmt, ... ) {
  int n_char_written;
  va_list args;
  va_start(args, fmt);
  n_char_written = vsprintf(ptr, fmt, args);
  va_end(args);
  ptr += n_char_written;
}
  
// CharBuffer::WriteDouble()
/// Write a double to the buffer with given format.
void CharBuffer::WriteDouble(const char *format, double dval) {
  int n_char_written;
  //printf("DEBUG:\tWriteDouble[%s, %lf]\n",format, dval);
  n_char_written = sprintf(ptr, format, dval);
  ptr += n_char_written;
}

// CharBuffer::WriteInteger()
/// Write an integer to the buffer with given format.
void CharBuffer::WriteInteger(const char *format, int ival) {
  int n_char_written;
  //printf("DEBUG:\tWriteInteger[%s, %i]\n",format, ival);
  n_char_written = sprintf(ptr, format, ival);
  ptr += n_char_written;
}

// CharBuffer::WriteString()
/// Write a string to the buffer with given format.
void CharBuffer::WriteString(const char *format, const char *sval) {
  int n_char_written;
  //printf("DEBUG:\tWriteString[%s, %s]\n",format, sval);
  n_char_written = sprintf(ptr, format, sval);
  ptr += n_char_written;
}

// CharBuffer::WriteDoubleXYZ()
/// Write the given double array of size 3 to the buffer with given format.
void CharBuffer::WriteDoubleXYZ(const char *format, double *XYZ) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, XYZ[0], XYZ[1], XYZ[2]);
  ptr += n_char_written;
}

// CharBuffer::NewLine()
/// Add a newline character to buffer.
void CharBuffer::NewLine() {
  ptr[0]='\n';
  ptr++;
}

// CharBuffer::Space()
/// Add space to buffer
void CharBuffer::Space() {
  ptr[0]=' ';
  ptr++;
}

int CharBuffer::Read(void *str, size_t numbytes) {
  size_t bytes_to_read;
  // Figure out actual size to copy
  size_t bytes_remaining = end_buffer - ptr;
  if (bytes_remaining < 1)
    return -1;
  else if (numbytes > bytes_remaining)
    bytes_to_read = bytes_remaining;
  else
    bytes_to_read = numbytes;

  // Read from ptr into str
  memcpy(str, ptr, bytes_to_read);
  // Advance ptr
  ptr += bytes_to_read;
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
  while (ptr < end_buffer && currentNum < num) {
    str[currentNum++] = *ptr;
    if (*ptr=='\n' || *ptr=='\0') {++ptr; break;}
    ++ptr;
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
char *BufferToDouble(char *buffer, double *X, int N, int width) {
  char *ptr;
  char number[64]; // Should not have to handle numbers wider than this!
  int i,atom;

//  number = (char*) malloc( (width+1) * sizeof(char));
  number[width]='\0';
  ptr=buffer;
  for (atom=0; atom<N; atom++) {
    i=0;
    while (i<width) {
      if (*ptr=='\n' || *ptr=='\r') ptr++;
      if (*ptr=='*') return NULL;//{free(number); return 1;}
      number[i++]=*ptr;
      ptr++;
    }
    //fprintf(stdout,"DEBUG: %i: %s\n",atom/3,number);
    X[atom] = atof(number);
  }

//  free(number);
  return ptr;
}

// DoubleToBuffer()
/** Given an array of double, format, and character width corresponding
  * to format, write coords in array to buffer. 
  * Return the position in the buffer after write. 
  */
char *DoubleToBuffer(char *buffer, double *X, int N, const char *format, 
                                 int width, int numCols) {
  int coord;
  char *ptr;

  ptr=buffer;
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
}

// BoxToBuffer()
/** Given an array of double[6], format, and character width corresponding
  * to format, write box coords to buffer.
  * Return the position in the buffer after write.
  */
char *BoxToBuffer(char *buffer, double *box, int numBox, 
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

