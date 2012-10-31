#ifndef INC_BUFFEREDFILE_H
#define INC_BUFFEREDFILE_H
#include "CpptrajFile.h"
/// Used to buffer text files that will be read line-by-line and/or in chunks
class BufferedFile : public CpptrajFile {
  public:
    BufferedFile();
    ~BufferedFile();

    size_t SetupFrameBuffer(int, int, int, int);
    size_t ResizeBuffer(int);
    int SeekToFrame(size_t);
    int ReadFrame();
    int WriteFrame();
    void GetDoubleAtPosition(double&,size_t,size_t);
    void BufferBeginOffset();
    void BufferBegin();
    void BufferToDouble(double*,int);
    void DoubleToBuffer(const double*,int, const char*);

    std::string GetLineUnbuffered();

    size_t FrameSize() { return frameSize_; }
    const char* Buffer() { return buffer_; }
  private:
    static const size_t DEFAULT_BUFFERSIZE = 16384;
    static const size_t LINE_BUF_SIZE = 1024;
    char linebuffer_[LINE_BUF_SIZE];

    char* buffer_;
    char* bufferPosition_;
    size_t frameSize_;
    size_t offset_;
    int Nelts_;
    int Ncols_;
    size_t eltWidth_;
};
#endif
