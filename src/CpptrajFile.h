#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include <string>
#include "FileIO.h"
// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
class CpptrajFile {
  public:
    CpptrajFile();
    virtual ~CpptrajFile(); // Virtual since class is inherited
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);

    /// Prepare file for reading. 
    int SetupRead(const char*, int);
    int SetupRead(std::string const&, int);
    int OpenRead(std::string const&);
    /// Prepare file for writing.
    int SetupWrite(const char*, int);
    int SetupWrite(std::string const&, int);
    int OpenWrite(std::string const&);
    /// Prepare file for appending. 
    int SetupAppend(const char*, int);
    /// Open file.
    int OpenFile();
    /// Close file.
    void CloseFile();
    /// Printf using the files Write routine.
    void Printf(const char*, ...);
    /// Printf using the files Write routine for the given rank.
    void Rank_printf(int, const char *, ...);
    /// Return true if the file is open
    bool IsOpen();
    /// Return the file name
    const char *Name();
    std::string FullPathName();
    /// Return the file name without the path
    const char *BaseName();
    /// Return the file extension
    std::string Extension();
    /// Return true if the file contains carriage returns in addition to newlines
    bool IsDos();
    /// Return true if the file is compressed.
    bool IsCompressed();

    int Gets(char* buf, int num) { return IO->Gets(buf, num); }
    int Write( char* buf, size_t num) { return IO->Write(buf, num); }
  protected:
    enum CompressType { NO_COMPRESSION, GZIP, BZIP2, ZIP };
    enum AccessType   { READ, WRITE, APPEND };

    FileIO* IO;                 ///< The interface to basic IO operations.
    AccessType access_;         ///< Access (Read, write, append)
    int isDos_;                 ///< 1 if CR present, need to count them as newlines
    off_t uncompressed_size_;   ///< If compressed, uncompressed file size
    off_t file_size_;           ///< Actual file size
    CompressType compressType_; ///< Type of compression
    std::string FileName_;      ///< Passed in filename
    int debug_;                 ///< Debug level
    bool isOpen_;               ///< If true, file is open and ready for IO.

  private:
    enum FileType {
      UNKNOWN_TYPE, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE
    };
    static const char FileTypeName[][13];
    static const char AccessTypeName[][2];
    static const size_t BUF_SIZE;

    FileType fileType_;         ///< File type (determines IO)
    std::string basefilename_;  ///< Filename minus any path
    // TODO: Make this a string
    char *filename_;            ///< Avoids constant calls to FileName_.c_str()
    std::string Ext_;           ///< Filename extension
    char printf_buffer_[1024];  ///< Used in Printf functions

    void Reset();
    void SetBaseFilename(const char*);
    int SetupFileIO();
    int ID_Type();
};
#endif
