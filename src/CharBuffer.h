#ifndef INC_CHARBUFFER_H
#define INC_CHARBUFFER_H
#include <cstddef> // for size_t
/*! \file CharBuffer.h
    \brief Contains a class and routines used for manipulating character buffers.
  */
// Class: CharBuffer
/// Used to manipulate character buffers.
class CharBuffer {
    char *buffer;      ///< The character buffer
    char *ptr;         ///< Current position in the buffer
    char *end_buffer;  ///< The end of the buffer
    size_t bufferSize; ///< Total size of the buffer
  public:
    CharBuffer();
    ~CharBuffer();

    char *Buffer() { return buffer; }

    void Allocate(size_t);
    void IncreaseSize(size_t);
    void Rewind();
    void Sprintf(const char *, ... ); 
    void WriteDouble(const char*,double);
    void WriteInteger(const char*,int);
    void WriteString(const char*,const char*);
    void WriteDoubleXYZ(const char*,double*);
    void NewLine();
    void Space();
    int Read(void *, size_t);
    int Gets(char *, int);
    size_t CurrentSize();
    /// Return allocated size of buffer
    size_t BufferSize() { return bufferSize; }
    /// Return a pointer to current position in buffer
    char *BufferPtr() { return ptr; }
};
// =============================================================================
char *BufferToDouble(char*,double*,int,int);
char *DoubleToBuffer(char*,double*,int,const char*,int,int);
char *BoxToBuffer(char*,double*,int,const char*,int);
#endif
