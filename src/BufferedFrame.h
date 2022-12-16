#ifndef INC_BUFFEREDFRAME_H
#define INC_BUFFEREDFRAME_H
#include <cstddef> // size_t
#include "CpptrajFile.h"
#include "TextFormat.h"
/// Used to buffer text files that will be read/written in formatted 'frames'.
/** Assumes the file is laid out as follows:
  * [<file header bytes>]
  * [<frame0 header bytes>]
  * <frame0>
  * [<frame1 header bytes>]
  * <frame1>
  * ...
  * Where <file header bytes> are optional non-frame bytes at the beginning
  * of the file, <frameX header bytes> are optional frame bytes at the
  * beginning of each frame, and <frameX> is the acual formatted frame
  * data.
  * NOTE: The preprocessor define CPPTRAJ_DEBUG_BUFFEREDFRAME can be used
  *       to provide more details on buffer write errors or overflows.
  */
class BufferedFrame : public CpptrajFile {
  public:
    /// CONSTRUCTOR
    BufferedFrame();
    /// DESTRUCTOR
    ~BufferedFrame();
    /// Set up buffer for read/write with # elts, elt width, elts per line
    size_t SetupFrameBuffer(int, int, int);
    /// Set up buffer for read/write. # elts, elt width, elts per line, additional bytes, offset.
    size_t SetupFrameBuffer(int, int, int, size_t, int);
    /// Set up buffer (primarily for writes). # elements, text format, elts per line.
    size_t SetupFrameBuffer(int, TextFormat const&, int);
    /// Change size of buffer by given delta
    size_t ResizeBuffer(int);

    /// Seek to a specific frame in the file
    int SeekToFrame(size_t);
    /// Set buffer position to beginning of buffer.
    void BufferBegin();
    /// Set buffer position to specified position.
    void BufferBeginAt(size_t);

    /// Attempt to read frameSize_ bytes.
    int AttemptReadFrame();
    /// Read frameSize_ bytes into buffer.
    bool ReadFrame();
    /// Convert contents of buffer to an array of double precision numbers
    void BufferToDouble(double*,int);
    /// \return Pointer to next element in buffer (null-terminated).
    const char* NextElement();

    /// Write frameSize_ bytes from buffer. // TODO should this be hidden in favor of FlushBuffer
    int WriteFrame();
    /// Convert array of double precision numbers to text in buffer.
    void DoubleToBuffer(const double*,int, const char*);
    /// Place integer into buffer and advance
    void IntToBuffer(int);
    /// Place double into buffer and advance
    void DblToBuffer(double);
    /// Place character string into buffer and advance
    void CharToBuffer(const char*);
    /// Write contents of buffer to file.
    void FlushBuffer();

    /// \return Current frame size
    size_t FrameSize()   const { return frameSize_; }
    /// \return Const pointer to buffer
    const char* Buffer() const { return buffer_;    }
    /// \return Total output file size for given number of frames.
    size_t OutputFileSize(unsigned int n) const { return offset_ + (frameSize_ * n); }
  private:
    size_t CalcFrameSize(int) const;
    inline void AdvanceCol();

    char* buffer_;         ///< Character buffer.
    char* bufferPosition_; ///< Current position in buffer_.
    size_t frameSize_;     ///< Total size of frame to read/write.
    size_t offset_;        ///< Non-frame bytes at the beginning of the file. Used in seeking.
    size_t memSize_;       ///< Current size of the buffer.
    size_t maxSize_;       ///< Actual size of the buffer in memory.
    int Ncols_;            ///< Number of columns in each frame.
    int col_;              ///< Current column (during writes).
    size_t eltWidth_;      ///< Width in chars of each text element in a frame.
    char saveChar_;        ///< For NextElement(); saved character at bufferPosition_.
    TextFormat writeFmt_;  ///< Write format to use for {Int|Dbl|Char}ToBuffer.
    int errorCount_;       ///< Number of errors since last write.
    int overflowCount_;    ///< Number of character overflows since last write.
};
#endif
