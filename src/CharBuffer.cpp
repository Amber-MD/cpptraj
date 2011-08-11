// CharBuffer
#include "CharBuffer.h"
#include <cstdio> // sprintf
#include <cstdlib> // atof

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

