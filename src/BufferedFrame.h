#ifndef INC_BUFFEREDFRAME_H
#define INC_BUFFEREDFRAME_H
#include "CpptrajFile.h"
#include "TextFormat.h"
/// Used to buffer text files that will be read in chunks
class BufferedFrame : public CpptrajFile {
  public:
    BufferedFrame();
    ~BufferedFrame();
    /// Set up buffer for read/write with # elts, elt width, elts per line
    size_t SetupFrameBuffer(int, int, int);
    /// Set up buffer for read/write. Can specify additional bytes in frame/offset
    size_t SetupFrameBuffer(int, int, int, size_t, int);
    /// Set up buffer (primarily for writes)
    size_t SetupFrameBuffer(int, TextFormat const&, int);
    size_t ResizeBuffer(int);
    int SeekToFrame(size_t);
    /// Attempt to read frameSize_ bytes.
    int AttemptReadFrame();
    /// Read frameSize_ bytes.
    bool ReadFrame();
    int WriteFrame();
    void GetDoubleAtPosition(double&,size_t,size_t);
    void BufferBegin();
    void BufferBeginAt(size_t);
    void AdvanceBuffer(size_t);
    void BufferToDouble(double*,int);
    void DoubleToBuffer(const double*,int, const char*);
    const char* NextElement();

    void IntToBuffer(int);
    void DblToBuffer(double);
    void CharToBuffer(const char*);
    void FlushBuffer();

    size_t FrameSize()   const { return frameSize_; }
    const char* Buffer() const { return buffer_;    }
    /// \return Total output file size for given number of frames.
    size_t OutputFileSize(unsigned int n) const { return offset_ + (frameSize_ * n); }
  private:
    size_t CalcFrameSize(int) const;
    inline void AdvanceCol();

    char* buffer_;         ///< Character buffer.
    char* bufferPosition_; ///< Position in buffer.
    size_t frameSize_;     ///< Total size of frame to read.
    size_t offset_;        ///< User specified offset, used in seeking.
    size_t memSize_;       ///< Total size of the buffer in memory.
    size_t maxSize_;       ///< Max size of the buffer in memory.
    int Ncols_;            ///< Number of columns, use to convert array to buffer.
    int col_;              ///< Current column (writes)
    size_t eltWidth_;      ///< Width of each element in the frame.
    char saveChar_;        ///< For NextElement(); saved character at bufferPosition_.
    TextFormat writeFmt_;  ///< Format to use for Int|Dbl|CharToBuffer
};
#endif
