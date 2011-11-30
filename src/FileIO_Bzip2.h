#ifndef INC_FILEIO_BZIP2_H
#define INC_FILEIO_BZIP2_H
#ifdef HASBZ2
// NOTE: bzlib.h has stdio.h. Does it matter that its not cstdio? 
#include "bzlib.h"
#include "FileIO.h"
// Class: FileIO_Bzip2
/// Bzip2 file IO.
class FileIO_Bzip2 : public FileIO {
    bool isBzread;
    FILE *fp;
    BZFILE *infile;
    int err;
    char *bzfilename;
    char *bzmode;
    off_t position;

    const char *BZerror();
  public:
    FileIO_Bzip2(); 
    ~FileIO_Bzip2(); 
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
