#ifndef INC_BZIP2FILE_H
#define INC_BZIP2FILE_H
// Bzip2 file IO
#ifdef HASBZ2
#include "BaseFileIO.h"
// NOTE: bzlib.h has stdio.h. Does it matter that its not cstdio? 
#include "bzlib.h" 
class Bzip2File : public BaseFileIO {
    bool isBzread;
    FILE *fp;
    BZFILE *infile;
    int err;
    char *bzfilename;
    char *bzmode;
    off_t position;

    const char *BZerror();
  public:
    Bzip2File(); 
    ~Bzip2File(); 
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
