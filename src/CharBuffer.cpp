// CharBuffer
#include "CharBuffer.h"
#include <cstdio> // sprintf
#include <cstdlib> // atof

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

/* CharBuffer::Allocate()
 */
void CharBuffer::Allocate(size_t sizeIn) {
  if (buffer!=NULL) delete[] buffer;
  bufferSize = sizeIn;
  buffer = new char[ bufferSize ];
  ptr = buffer;
}

/* CharBuffer::CurrentSize()
 */
size_t CharBuffer::CurrentSize() {
  return (size_t) (ptr - buffer);
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

/* CharBuffer::WriteDouble()
 */
void CharBuffer::WriteDouble(const char *format, double dval) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, dval);
  ptr += n_char_written;
}

/* CharBuffer::WriteInteger()
 */
void CharBuffer::WriteInteger(const char *format, int ival) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, ival);
  ptr += n_char_written;
}

/* CharBuffer::WriteString()
 */
void CharBuffer::WriteString(const char *format, const char *sval) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, sval);
  ptr += n_char_written;
}

/* CharBuffer::WriteString()
 */
void CharBuffer::WriteString(const char *sval) {
  int n_char_written;
  n_char_written = sprintf(ptr, "%s", sval);
  ptr += n_char_written;
}


/* CharBuffer::WriteStringN()
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

/* CharBuffer::WriteDoubleXYZ()
 */
void CharBuffer::WriteDoubleXYZ(const char *format, double *XYZ) {
  int n_char_written;
  n_char_written = sprintf(ptr, format, XYZ[0], XYZ[1], XYZ[2]);
  ptr += n_char_written;
}

void CharBuffer::NewLine() {
  ptr[0]='\n';
  ptr++;
}

// =============================================================================
/* BufferToDouble()
 * Store character buffer containing XYZ coords with format 
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

/* DoubleToBuffer()
 * Given an array of double, format, and character width corresponding
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

/* BoxToBuffer()
 * Given an array of double[6], format, and character width corresponding
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

