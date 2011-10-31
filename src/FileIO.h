#ifndef INC_FILEIO_H
#define INC_FILEIO_H
#include <sys/types.h> // For off_t
/// Class: FileIO 
/// Base class for handling basic file IO. File types will inherit this
/// class and replace basic IO functions (open, close, read etc) with
/// their own.
/// NOTES:
///   File pointers for each file type will be declared in their own class.
///   Use off_t for seeking for better compatibility with large files.
///   Thus far it has only been necessary to seek using SEEK_SET, so that
///   argument will be removed from Seek() for now. 
class FileIO {
  public:
    virtual ~FileIO() {}
    virtual int Open(const char *, const char *) { return 0;  }
    virtual int Close()                          { return 0;  }
    virtual int Read(void *, size_t, size_t)     { return 0;  }
    virtual int Write(void *, size_t, size_t)    { return 0;  }
    virtual int Seek(off_t)                      { return 0;  }
    virtual int Rewind()                         { return 0;  }
    virtual off_t Tell()                         { return 0;  }
    virtual int Gets(char *, int)                { return 0;  }
    // Only needed for compressed files
    virtual off_t Size(char *)                   { return 0;  }
    // Only needed for MPI file
    virtual int SetSize(long int)                { return 0;  }
    // Printf routines that call Write directly
    void Printf(const char*, ...);
    void Rank_printf(int, const char *, ...);
};
#endif
