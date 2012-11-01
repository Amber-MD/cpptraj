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
    void BufferBeginAt(size_t);
    void BufferToDouble(double*,int);
    void DoubleToBuffer(const double*,int, const char*);

    std::string GetLineUnbuffered();

    size_t FrameSize() { return frameSize_; }
    const char* Buffer() { return buffer_; }
  private:
    static const size_t DEFAULT_BUFFERSIZE = 16384;

    char* buffer_;
    char* bufferPosition_;
    size_t frameSize_;     ///< Total size of frame to read.
    size_t offset_;        ///< User specified offset, used in seeking.
    int Ncols_;            ///< Number of columns, use to convert array to buffer.
    size_t eltWidth_;      ///< Width of each element in the frame.
};
#endif
