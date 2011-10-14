#ifndef INC_CHARBUFFER_H
#define INC_CHARBUFFER_H
/// CharBuffer
/// Used to manipulate character buffers.
/*
class CharBuffer {
    char *buffer;      // The character buffer
    char *ptr;         // Current position in the buffer
    size_t bufferSize; // Total size of the buffer
  public
    CharBuffer();
    ~CharBuffer();

    int Setup(size_t);
    int Sprintf(const char *, ...);
};
*/
char *BufferToDouble(char*,double*,int,int);
char *DoubleToBuffer(char*,double*,int,const char*,int,int);
char *BoxToBuffer(char*,double*,int,const char*,int);
#endif
