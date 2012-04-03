#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include <string>
#include "FileIO.h"
#ifdef USE_CHARBUFFER
#include "CharBuffer.h"
#endif
/* Experimental Compiler Defines:
 * -USE_CHARBUFFER: Use CharBuffer to buffer an entire file.
 */ 
// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
class CpptrajFile {
  public:
    /** This is the interface to basic IO operations. It is exposed in 
      * order to speed up operations that are frequently called, like
      * Reads/Writes.
      */
    FileIO *IO; 

    CpptrajFile();
    virtual ~CpptrajFile(); // Virtual since class is inherited
    CpptrajFile(const CpptrajFile&);
    CpptrajFile &operator=(const CpptrajFile &);

    /// Prepare file for reading. 
    int SetupRead(char*, int);
    /// Prepare file for writing.
    int SetupWrite(char*, int);
    /// Prepare file for appending. 
    int SetupAppend(char*, int);
    /// Open file.
    int OpenFile();
    /// Close file.
    void CloseFile();
    /// Printf using the files Write routine.
    void Printf(const char*, ...);
    /// Printf using the files Write routine for the given rank.
    void Rank_printf(int, const char *, ...);
    // frameBuffer routines
    void BufferBegin(size_t);
    void BufferBegin();
    void BufferToDouble(double *, int, int);
    void DoubleToBuffer(double *, int, const char*, int, int);
    void BoxToBuffer(double *, int, const char*, int);
#   ifdef USE_CHARBUFFER
    int OpenFileBuffered();
    //int ReadBuffered();
    int Gets(char*, int);
    void Rewind();
    int Read(void*,size_t);
#   endif
    /// Return true if the file is open
    bool IsOpen();
    /// Return the file name
    const char *Name();
    std::string FullPathName();
    /// Return the file name without the path
    const char *BaseName();
    /// Return the file extension
    std::string Extension();
    /// Return true if filename or basefilename matches
    bool FilenameIs(const std::string&);
    bool FilenameIs(char*);
    /// Return true if the file contains carriage returns in addition to newlines
    bool IsDos();
    /// Return true if the file is compressed.
    bool IsCompressed();
  protected:
    enum CompressType {
      NO_COMPRESSION, GZIP, BZIP2, ZIP
    };
    enum AccessType {
      READ, WRITE, APPEND
    };
    AccessType access_;         ///< Access (Read, write, append)
    int isDos_;                 ///< 1 if CR present, need to count them as newlines
    off_t uncompressed_size_;   ///< If compressed, uncompressed file size
    off_t file_size_;           ///< Actual file size
    CompressType compressType_; ///< Type of compression
    std::string FileName_;      ///< Passed in filename
    int debug_;                 ///< Debug level
    bool isOpen_;               ///< If true, file is open and ready for IO.
    // frameBuffer vars
    // TODO: Make private
    char *frameBuffer_;         ///< Used to buffer frames for Amber traj/restart
    char *bufferPosition_;      ///< Position in frameBuffer
    size_t frameSize_;          ///< Size of frameBuffer
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

#   ifdef USE_CHARBUFFER
    CharBuffer c_buffer_;
#   endif

    void Reset();
    void SetBaseFilename(char*);
    int SetupFileIO();
    int ID_Type();
};
#endif
