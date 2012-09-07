#include "FileBuffer.h"
#include "CpptrajStdio.h"

FileBuffer::FileBuffer() :
  IO_(0),
  endlinebuffer_(linebuffer_ + LINE_BUF_SIZE),
  total_read_(0),
  readbuffer_(0),
  lineptr_(linebuffer_),
  ptr_(0),
  endbuffer_(0) 
{}

FileBuffer::FileBuffer(FileIO* IOin, int sizeIn) :
  IO_(IOin),
  endlinebuffer_(linebuffer_ + LINE_BUF_SIZE),
  total_read_(0),
  readbuffer_(0),
  lineptr_(linebuffer_),
  ptr_(0),
  endbuffer_(0),
  progress_(sizeIn)
{
  readbuffer_ = new char[ DEFAULT_CHUNKSIZE ];
  ptr_ = readbuffer_;
  endbuffer_ = ptr_; // This guarantees on first read buffer will be filled
}

FileBuffer::~FileBuffer() {
  if (readbuffer_!=0) delete[] readbuffer_;
}

const char* FileBuffer::NextLine() {
  // Reset line buffer
  lineptr_ = linebuffer_;
  // Get next line from chunk.
  while ( lineptr_ < endlinebuffer_ ) {
    // Fill buffer if needed.
    if (endbuffer_ == ptr_) {
      int Nread = IO_->Read(readbuffer_, DEFAULT_CHUNKSIZE);
      if (Nread < 1) return 0;
      total_read_ += Nread;
      progress_.Update( total_read_ );
      ptr_ = readbuffer_;
      endbuffer_ = readbuffer_ + (size_t)Nread;
    }
    // Fill line buffer
    *(lineptr_++) = *ptr_;
    if (*ptr_ == '\n') {
      *lineptr_ = '\0';
      // Position ptr at next char
      ++ptr_;
      return linebuffer_;
    }
    ++ptr_;
  }
  mprinterr("Error: FileBuffer: blowing line buffer (> %zu bytes)\n", LINE_BUF_SIZE);
  linebuffer_[LINE_BUF_SIZE - 1] = '\0';
  return linebuffer_;
}
