#ifndef INC_FILEIO_STD_H
#define INC FILEIO_STD_H
#include "FileIO.h"
#include <cstdio> // For FILE
/// Class: FileIO_Std
/// Standard file IO
class FileIO_Std : public FileIO {
    FILE *fp;
    bool isStdout;
  public:
    FileIO_Std();
    ~FileIO_Std();
    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t, size_t );
    int Write(void *, size_t, size_t);  
    int Seek(off_t);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
};
#endif
