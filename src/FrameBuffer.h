#ifndef INC_FRAMEBUFFER_H
#define INC_FRAMEBUFFER_H
#include <cstddef> // size_t
class FrameBuffer {
  public : 
    FrameBuffer();
    ~FrameBuffer();
    FrameBuffer(int,int,int,bool);
    FrameBuffer(const FrameBuffer&);
    FrameBuffer& operator=(const FrameBuffer&);
    char& operator[](int idx) { return frameBuffer_[idx]; } // NOTE: No bounds check

    void BufferBegin(size_t);
    void BufferBegin();
    void BufferToDouble(double *, int, int);
    void DoubleToBuffer(const double *, int, const char*, int, int);
    void BoxToBuffer(const double *, int, const char*, int);

    char* Buffer()       { return frameBuffer_;    }
    char* CurrentPos()   { return bufferPosition_; }
    size_t CurrentSize() { return (bufferPosition_ - frameBuffer_); }
    size_t ReadSize()    { return frameSize_ - 1; } // No NULL char
    size_t FrameSize()   { return frameSize_;     }
    void MovePos(int offset) { bufferPosition_ += offset; }
  private :
    char *frameBuffer_;
    char *bufferPosition_;
    size_t frameSize_;
};
#endif
