#ifndef INC_FRAMEBUFFER_H
#define INC_FRAMEBUFFER_H
#include <cstddef> // size_t
class FrameBuffer {
  public : 
    FrameBuffer();
    ~FrameBuffer();
    FrameBuffer(const FrameBuffer&);
    FrameBuffer& operator=(const FrameBuffer&);

    void BufferBegin(size_t);
    void BufferBegin();
    void BufferToDouble(double *, int, int);
    void DoubleToBuffer(double *, int, const char*, int, int);
    void BoxToBuffer(double *, int, const char*, int);
  protected :
    char *frameBuffer_;
    char *bufferPosition_;
    size_t frameSize_;
};
#endif
