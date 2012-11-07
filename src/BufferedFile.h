#ifndef INC_BUFFEREDFILE_H
#define INC_BUFFEREDFILE_H
#include "CpptrajFile.h"
/// Used to buffer text files that will be read line-by-line and/or in chunks
class BufferedFile : public CpptrajFile {
  public:
    BufferedFile();
    ~BufferedFile();

    int SetupBuffer();
    const char* BufferedLine();
    int TokenizeLine(const char*);
    const char* NextToken();

    size_t SetupFrameBuffer(int, int, int);
    size_t SetupFrameBuffer(int, int, int, size_t, int);
    size_t ResizeBuffer(int);
    int SeekToFrame(size_t);
    int ReadFrame();
    int WriteFrame();
    void GetDoubleAtPosition(double&,size_t,size_t);
    void BufferBegin();
    void BufferBeginAt(size_t);
    void AdvanceBuffer(size_t);
    void BufferToDouble(double*,int);
    void DoubleToBuffer(const double*,int, const char*);

    size_t FrameSize() { return frameSize_; }
    const char* Buffer() { return buffer_; }
  private:
    static const size_t DEFAULT_BUFFERSIZE = 16384;

    size_t CalcFrameSize(int);

    char* buffer_;         ///< Character buffer
    char* bufferPosition_; ///< Position in buffer
    size_t frameSize_;     ///< Total size of frame to read.
    size_t offset_;        ///< User specified offset, used in seeking.
    int Ncols_;            ///< Number of columns, use to convert array to buffer.
    size_t eltWidth_;      ///< Width of each element in the frame.

    /// Array of pointers to beginning of tokens in linebuffer
    char* tokens_[BUF_SIZE];
    char** tokenptr_;      ///< Position in tokens
    char savechar_;        ///< Saved last char of current token
    char* lineptr_;        ///< Position in linebuffer
    char* endlinebuffer_;  ///< End position of linebuffer
    char* endbuffer_;      ///< End position of buffer
};
#endif
