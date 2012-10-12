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

    void BufferBegin(size_t);
    void BufferBegin();
    void BufferToDouble(double *, int, int);
    void DoubleToBuffer(const double *, int, const char*, int, int);
    void BoxToBuffer(const double *, int, const char*, int);

    char* Buffer() { return frameBuffer_; }
    size_t CurrentSize() { return (bufferPosition_ - frameBuffer_); }
  protected :
    char *frameBuffer_;
    char *bufferPosition_;
    size_t frameSize_;
};
#endif
