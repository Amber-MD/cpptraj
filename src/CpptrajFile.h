#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include "FileName.h" 
#include "FileIO.h"
#ifdef MPI
# include "Parallel.h"
#endif
/// Class to abstract handling of basic file routines.
class CpptrajFile {
  public:
    enum AccessType   { READ=0, WRITE, APPEND, UPDATE };
    enum CompressType { NO_COMPRESSION=0, GZIP, BZIP2, ZIP };
    enum FileType { UNKNOWN_TYPE=0, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE };

    CpptrajFile();
    virtual ~CpptrajFile(); // Virtual since class is inherited
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);
    /// Set debug level.
    void SetDebug(int d) { debug_ = d; }
    /// Open file with current access.
    int OpenFile() { return OpenFile( access_ ); }
    /// Open file with given access.
    int OpenFile(AccessType);
    /// Close file.
    void CloseFile();
    // -------------------------------------------
    /// Printf using the Write routine.
    void Printf(const char*, ...);
    /// Get next line as a string
    std::string GetLine();
    /// Get next line and return pointer to raw buffer
    const char* NextLine();
    // -------------------------------------------
    /// Set up and open file for reading.
    int OpenRead(FileName const&);
    /// Set up and open file for writing.
    int OpenWrite(FileName const&);
    /// Open file previously set up for write with given numeric suffix/extension prefix.
    int OpenWriteNumbered(int,bool);
    /// Set up and open file for appending.
    int OpenAppend(FileName const&);
    // -------------------------------------------
    /// Prepare file for reading. 
    int SetupRead(FileName const&, int);
    /// Prepare file for writing.
    int SetupWrite(FileName const&, int);
    /// Prepare file of given type for writing
    int SetupWrite(FileName const&, FileType, int);
    /// Prepare file for appending. 
    int SetupAppend(FileName const&, int);
    // -------------------------------------------
#   ifdef MPI
    int ParallelOpenFile(AccessType, Parallel::Comm const&);
    int ParallelOpenFile(Parallel::Comm const& c) { return ParallelOpenFile(access_, c); }
    // -------------------------------------------
#   endif
    /// \return the access file is currently set up for.
    AccessType Access()         const { return access_;               }
    /// \return the compression type
    CompressType Compression()  const { return compressType_;         }
    /// \return true if the file is open
    bool IsOpen()               const { return isOpen_;               }
    /// \return true if the file is a stream
    bool IsStream()             const { return isStream_;             }
    /// \return file name class.
    const FileName& Filename()  const { return fname_;                }
    /// \return 1 if the file contains carriage returns in addition to newlines
    int IsDos()                 const { return isDos_;                }
    /// \return File size
    off_t FileSize()            const { return file_size_;            }
    /// \return true if the file is compressed.
    bool IsCompressed()         const { return !(compressType_ == NO_COMPRESSION); }
    /// \return uncompressed file size (just size if file is not compressed).
    off_t UncompressedSize() const;
    int Gets(char* buf, int num)           { return IO_->Gets(buf, num);  }
    int Write(const void* buf, size_t num) { return IO_->Write(buf, num); }
    int Read(void* buf, size_t num)        { return IO_->Read(buf, num);  }
    int Seek(off_t offset)                 { return IO_->Seek(offset);    }
    int Rewind()                           { return IO_->Rewind();        }
    int Flush()                            { return IO_->Flush();         }
    off_t Tell()                           { return IO_->Tell();          }
    static const size_t BUF_SIZE = 1024;
  protected: // Protected for PDBfile/Mol2File
    char linebuffer_[BUF_SIZE]; ///< Used in Printf/GetLine functions
  private:
    static const char* FileTypeName[];
    static const char* AccessTypeName[];
    FileIO* IO_;                ///< The interface to basic IO operations.
    AccessType access_;         ///< Access (Read, write, append)
    int isDos_;                 ///< 1 if CR present, need to count them as newlines
    off_t uncompressed_size_;   ///< If compressed, uncompressed file size
    off_t file_size_;           ///< Actual file size
    CompressType compressType_; ///< Type of compression
    int debug_;                 ///< Debug level
    bool isOpen_;               ///< If true, file is open and ready for IO.
    bool isStream_;             ///< If true, file IO is to/from a stream.
    FileType fileType_;         ///< File type (determines IO)
    FileName fname_;            ///< Holds full and base file name + any extensions.

    void Reset();
    static FileIO* SetupFileIO(FileType);
    int ID_Type(const char*);
};
#endif
