// CharBuffer
#include <cstdio>  // sprintf
#include <cstdlib> // atof
#include <cstdarg> // va_list
#include <cstring> // memcpy
#include "CharBuffer.h"

// CONSTRUCTOR
CharBuffer::CharBuffer() {
  buffer = NULL;
  ptr = NULL;
  bufferSize = 0;
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
}

// CharBuffer::AddCharString()
/// Add NULL terminated char string to buffer, increasing size accordingly.
void CharBuffer::AddCharString(char *inputString) {
  // Calculate size necessary to accomodate input string
  size_t inputStringSize = strlen(inputString);
  size_t newBufferSize = inputStringSize + CurrentSize() + 1;
  // Reallocate if the resulting buffer size greater than current allocation
  if (newBufferSize > bufferSize) {
    char *tempbuffer = new char[ newBufferSize ];
    if (buffer!=NULL) {
      memcpy(tempbuffer, buffer, bufferSize * sizeof(char));
      delete[] buffer;
    }
    buffer = tempbuffer;
    ptr = buffer + bufferSize;
    bufferSize = newBufferSize;
  }
  // Copy inputstring to buffer.
  memcpy(ptr, inputString, inputStringSize * sizeof(char));
  // Add terminating NULL character
  ptr += inputStringSize;
  *ptr = '\0';
}

// CharBuffer::CurrentSize()
/// Return the size of the data that has been currently written to the buffer.
size_t CharBuffer::CurrentSize() {
  //printf("DEBUG:\tCurrentSize=%lu, Allocd for %lu\n",(size_t) (ptr - buffer),bufferSize);
  return (size_t) (ptr - buffer);
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
  
/* CharBuffer::WriteStringBuffer()
 */
/*
void CharBuffer::WriteStringBuffer(char *sval) {
  size_t len = 0;
  char *s_ptr = sval;
  while (*s_ptr!='\0') {
    s_ptr++;
    len++;
  }
  size_t currentSize = CurrentSize();
  size_t newSize = currentSize + len + 1;
  // If not enough room to add sval to buffer, reallocate buffer
  if ( newSize > bufferSize ) {
    char *newBuffer = new char[ newSize ];
    memcpy(newBuffer, buffer, currentSize * sizeof(char));
    delete[] buffer;
    buffer = newBuffer;
    // Reposition ptr to same place it was in old buffer
    ptr = buffer + currentSize;
  }
  strcpy(ptr, sval);
  ptr += (len + 1); // +1 for NULL
}*/

// CharBuffer::WriteDouble()
/// Write a double to the buffer with given format.
void CharBuffer::WriteDouble(const char *format, double dval) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, dval);
  ptr += n_char_written;
}

// CharBuffer::WriteInteger()
/// Write an integer to the buffer with given format.
void CharBuffer::WriteInteger(const char *format, int ival) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, ival);
  ptr += n_char_written;
}

// CharBuffer::WriteString()
/// Write a string to the buffer with given format.
void CharBuffer::WriteString(const char *format, const char *sval) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, sval);
  ptr += n_char_written;
}

/* CharBuffer::WriteString()
 */
/*void CharBuffer::WriteString(const char *sval) {
  int n_char_written;
  n_char_written = sprintf(ptr, "%s", sval);
  ptr += n_char_written;
}*/

// CharBuffer::WriteStringN()
/** Write a string to the buffer with given width. If the string is larger
  * than <width>, truncate it. If leftAlign is specified align the string
  * to the left, preceded by '#'. Also write a space, so that the total
  * number of characters written is always <width> + 1.
  */
void CharBuffer::WriteStringN(char *sval, int width, bool leftAlign) {
  int c, n_char_written, actualWidth;
  actualWidth = width;
  if (leftAlign) actualWidth--;
  char *temps = new char[ actualWidth + 1];
  for (c = 0; c < actualWidth; c++) { 
    temps[c] = sval[c];
    if (temps[c]=='\0') break;
  }
  temps[c]='\0';
  if (leftAlign)
    n_char_written = sprintf(ptr, "#%-*s ",actualWidth,temps);
  else
    n_char_written = sprintf(ptr, " %*s",actualWidth,temps);
  //printf("DEBUG:\tWriteStringN(%s) width=%i n_char_written=%i\n",sval,width,n_char_written);
  delete[] temps;
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

