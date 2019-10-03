#ifndef INC_FILEIO_BZIP2_H
#define INC_FILEIO_BZIP2_H
#ifdef HASBZ2
// NOTE: bzlib.h has stdio.h. Does it matter that its not cstdio? 
#include <bzlib.h>
#include "FileIO.h"
// Class: FileIO_Bzip2
/// Bzip2 file IO.
class FileIO_Bzip2 : public FileIO {
  public:
    FileIO_Bzip2(); 
    ~FileIO_Bzip2();
    int OpenStream(StreamType) { return 1; } 
    int Open(const char *, const char *);    
    int Close();
    off_t Size(const char *);
    int Read(void *, size_t );
    int Write(const void *, size_t);
    int Flush() { return 0; } // bzflush is not part of the bzip standard
    int Seek(off_t);
    int Rewind();  
    off_t Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
    int SetSize(long int) { return 0; }
  private:
    static const char *BZerror(int);

    FILE *fp_;         ///< The system file handle.
    BZFILE *infile_;   ///< The bzip file handle.
    char *bzfilename_; ///< The file name, saved in case we need to rewind.
    char *bzmode_;     ///< The file mode, saved in case we need to rewind.
    off_t position_;   ///< The current position in the file, for Tell().
    int err_;          ///< The current error status.
    bool eofStat_;     ///< True if BZ_STREAM_END encountered during read.
    bool isBzread_;    ///< True if reading, false otherwise.
};
#endif
#endif
