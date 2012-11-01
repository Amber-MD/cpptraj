#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include "FileName.h" 
#include "FileIO.h"
// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
class CpptrajFile {
  public:
    enum AccessType   { READ, WRITE, APPEND };
    enum CompressType { NO_COMPRESSION, GZIP, BZIP2, ZIP };
    enum FileType { UNKNOWN_TYPE, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };

    CpptrajFile();
    virtual ~CpptrajFile(); // Virtual since class is inherited
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);
    /// Set up and open file for reading.
    int OpenRead(std::string const&);
    /// Prepare file for reading. 
    int SetupRead(std::string const&, int);
    /// Set up and open file for writing
    int OpenWrite(std::string const&);
    /// Prepare file for writing.
    int SetupWrite(std::string const&, int);
    /// Prepare file of given type for writing
    int SetupWrite(std::string const&, FileType, int);
    /// Set up and open file for appending.
    int OpenAppend(std::string const&);
    /// Prepare file for appending. 
    int SetupAppend(std::string const&, int);
    /// Switch file access mode
    int SwitchAccess(AccessType);
    /// Open file.
    int OpenFile();
    /// Close file.
    void CloseFile();
    /// Printf using the files Write routine.
    void Printf(const char*, ...);
    /// Printf using the files Write routine for the given rank.
    void Rank_printf(int, const char *, ...);
    /// Return the access file is currently set up for.
    AccessType Access()               { return access_;               }
    /// Return the compression type
    CompressType Compression()        { return compressType_;         }
    /// Return true if the file is open
    bool IsOpen()                     { return isOpen_;               }
    /// Return the file name with full path.
    const char* FullFileStr()         { return fname_.Full().c_str(); }
    /// String version of file name with full path.
    std::string const& FullFileName() { return fname_.Full();         }
    /// Return the file name without the path.
    const char* BaseFileStr()         { return fname_.Base().c_str(); }
    /// String version of file name without the path.
    std::string const& BaseFileName() { return fname_.Base();         }
    /// Return the file extension
    std::string const& Extension()    { return fname_.Ext();          }
    /// Return 1 if the file contains carriage returns in addition to newlines
    int IsDos() { return isDos_; }
    /// Return true if the file is compressed.
    bool IsCompressed();
    /// Return uncompressed file size (just size if file is not compressed).
    off_t UncompressedSize();
    FileIO* IOptr()                  { return IO_; } // TODO: Remove 
    int Gets(char* buf, int num)     { return IO_->Gets(buf, num);  }
    int Write(void* buf, size_t num) { return IO_->Write(buf, num); }
    int Read(void* buf, size_t num)  { return IO_->Read(buf, num);  }
    int Seek(off_t offset)           { return IO_->Seek(offset);    }
    int Rewind()                     { return IO_->Rewind();        }
  protected:
    static const size_t BUF_SIZE = 1024;
    char linebuffer_[BUF_SIZE]; ///< Used in Printf functions
  private:
    static const char FileTypeName[][13];
    static const char AccessTypeName[][2];

    FileIO* IO_;                ///< The interface to basic IO operations.
    AccessType access_;         ///< Access (Read, write, append)
    int isDos_;                 ///< 1 if CR present, need to count them as newlines
    off_t uncompressed_size_;   ///< If compressed, uncompressed file size
    off_t file_size_;           ///< Actual file size
    CompressType compressType_; ///< Type of compression
    int debug_;                 ///< Debug level
    bool isOpen_;               ///< If true, file is open and ready for IO.
    FileType fileType_;         ///< File type (determines IO)
    FileName fname_;            ///< Holds full and base file name + any extensions.

    void Reset();
    static FileIO* SetupFileIO(FileType);
    int ID_Type(const char*);
};
#endif
