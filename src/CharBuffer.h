#ifndef INC_CHARBUFFER_H
#define INC_CHARBUFFER_H
#include <cstddef> // for size_t
/*! \file CharBuffer.h
    \brief Contains a class and routines used for manipulating character buffers.
  */
// Class: CharBuffer
/// Used to manipulate character buffers.
class CharBuffer {
  public:
    CharBuffer();
    ~CharBuffer();
    CharBuffer(const CharBuffer&);
    CharBuffer &operator=(const CharBuffer&);
    /// Return a pointer to the beginning of the buffer
    char *Buffer();
    /// Return NULL terminated string
    const char *c_str();
    /// Return allocated size of buffer
    size_t BufferSize();
    /// Return a pointer to current position in buffer
    char *BufferPtr();
    /// Return size of data that is currently written to the buffer.
    size_t CurrentSize();
    /// Allocate space for buffer. Wipes out existing contents.
    void Allocate(size_t);
    /// Increase size of the buffer. Keeps current contents, updates ptr.
    void IncreaseSize(size_t);
    /// Increase size of the buffer. Keeps current contents.
    void Reallocate(size_t);
    /// Set pointer to beginning of buffer.
    void Rewind();
    /// Write data to the buffer with printf-like syntax.
    void Sprintf(const char *, ... ); 

    void DumpBuffer();
    void WriteDouble(const char*,double);
    void WriteInteger(const char*,int);
    void WriteString(const char*,const char*);
    void WriteXY(const char *, double, double);
    void WriteDoubleXYZ(const char*,double*);
    void NewLine();
    void Space();
    int Read(void *, size_t);
    int Gets(char *, int);
  private:
    char *buffer_;      ///< The character buffer
    char *ptr_;         ///< Current position in the buffer
    char *end_buffer_;  ///< The end of the buffer
    size_t bufferSize_; ///< Total size of the buffer
};
// =============================================================================
//char *BufferToDouble(char*,double*,int,int);
//char *DoubleToBuffer(char*,double*,int,const char*,int,int);
//char *BoxToBuffer(char*,double*,int,const char*,int);
#endif
