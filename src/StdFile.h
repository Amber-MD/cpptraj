#ifndef INC_STDFILE_H
#define INC STDFILE_H
// Standard file IO
#include "BaseFileIO.h" // cstdio
class StdFile : public BaseFileIO {
    FILE *fp;
    bool isStdout;
  public:
    StdFile();
    ~StdFile();
    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t, size_t );
    int Write(void *, size_t, size_t);  
    int Seek(off_t, int);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
};
#endif
