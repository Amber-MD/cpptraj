#ifndef INC_FILEBUFFER_H
#define INC_FILEBUFFER_H
#include "FileIO.h"
#include "ProgressBar.h"
class FileBuffer {
  public:
    FileBuffer();
    FileBuffer(FileIO*,int);
    ~FileBuffer();

    const char* NextLine();
  private:
    static const size_t DEFAULT_CHUNKSIZE = 16384;
    static const size_t LINE_BUF_SIZE = 1024;
    char linebuffer_[LINE_BUF_SIZE];
    FileIO* IO_;
    const char* endlinebuffer_;
    int total_read_;
    char* readbuffer_;
    char* lineptr_;
    char* ptr_;
    char* endbuffer_;
    ProgressBar progress_;
};
#endif 
