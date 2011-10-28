#ifndef INC_CHARBUFFER_H
#define INC_CHARBUFFER_H
/// Class: CharBuffer
/// Used to manipulate character buffers.
#include <cstddef>
class CharBuffer {
    char *buffer;      // The character buffer
    char *ptr;         // Current position in the buffer
    size_t bufferSize; // Total size of the buffer
  public:
    CharBuffer();
    ~CharBuffer();

    char *Buffer() { return buffer; }

    void Allocate(size_t);
    void WriteDouble(const char*,double);
    void WriteInteger(const char*,int);
    void WriteString(const char*,const char*);
    void WriteStringN(char *, int, bool);
    void WriteDoubleXYZ(const char*,double*);
    void NewLine();
    size_t CurrentSize();
};
// =============================================================================
char *BufferToDouble(char*,double*,int,int);
char *DoubleToBuffer(char*,double*,int,const char*,int,int);
char *BoxToBuffer(char*,double*,int,const char*,int);
#endif
