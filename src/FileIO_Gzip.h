#ifndef INC_FILEIO_GZIP_H
#define INC_FILEIO_GZIP_H
#ifdef HASGZ
#include "zlib.h"
#include "FileIO.h" 
// Class: FileIO_Gzip
/// Gzip file IO
class FileIO_Gzip : public FileIO {
    gzFile fp;
  public:
    FileIO_Gzip(); 
    ~FileIO_Gzip(); 
    int Open(const char *, const char *);    
    int Close();
    off_t Size(char *);
    int Read(void *, size_t, size_t );
    int Write(void *, size_t, size_t);  
    int Seek(off_t);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
};
#endif
#endif
