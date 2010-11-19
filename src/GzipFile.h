#ifndef INC_GZIPFILE_H
#define INC_GZIPFILE_H
#ifdef HASGZ
// Gzip file IO
#include "BaseFileIO.h" // cstdio
#include "zlib.h"
class GzipFile : public BaseFileIO {
    gzFile fp;
  public:
    GzipFile(); 
    ~GzipFile(); 
    int Open(const char *, const char *);    
    int Close();
    long long int Size(char *);
    int Read(void *, size_t, size_t );
    int Write(void *, size_t, size_t);  
    int Seek(off_t, int);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
};
#endif
#endif
