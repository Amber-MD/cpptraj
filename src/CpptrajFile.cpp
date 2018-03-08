// CpptrajFile
#include <sys/stat.h> // stat
#include <cstdio>     // vsprintf
#include <cstring>    // strlen, strerror 
#include <cerrno>     // errno
#include <cstdarg>    // va_X functions
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // AppendNumber
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "FileIO_Mpi.h"
#  include "FileIO_MpiShared.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

const char* CpptrajFile::FileTypeName[] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE", "MPISHARED"
};

const char* CpptrajFile::AccessTypeName[] = {
  "read", "write", "append", "update" };

// CONSTRUCTOR
CpptrajFile::CpptrajFile() :
  IO_(0),
  access_(READ),
  isDos_(0),
  uncompressed_size_(0UL),
  file_size_(0UL),
  compressType_(NO_COMPRESSION),
  debug_(0),
  isOpen_(false),
  isStream_(false),
  fileType_(STANDARD)
{ }

// Copy Constructor
CpptrajFile::CpptrajFile(const CpptrajFile &rhs) :
  IO_(0),
  access_(rhs.access_),
  isDos_(rhs.isDos_),
  uncompressed_size_(rhs.uncompressed_size_),
  file_size_(rhs.file_size_),
  compressType_(rhs.compressType_),
  debug_(rhs.debug_),
  // Even if file is open, copy it closed?
  isOpen_(false),
  isStream_(rhs.isStream_),
  fileType_(rhs.fileType_),
  fname_(rhs.fname_)
{
  // Set up the IO object
  // NOTE: Should probably throw an exception if this fails.
  if (rhs.IO_ != 0) 
    IO_ = SetupFileIO( fileType_ );
}

// Assignment
CpptrajFile &CpptrajFile::operator=(const CpptrajFile &rhs) {
  if (this != &rhs) { 
    // Deallocate
    CloseFile();
    if (IO_ != 0) delete IO_;
    // Allocate and copy
    debug_ = rhs.debug_;
    isStream_ = rhs.isStream_;
    access_ = rhs.access_;
    uncompressed_size_ = rhs.uncompressed_size_;
    file_size_ = rhs.file_size_;
    fileType_ = rhs.fileType_;
    fname_ = rhs.fname_;
    compressType_ = rhs.compressType_;
    isDos_ = rhs.isDos_;
    // Set up the IO object
    // NOTE: Should probably throw an exception if this fails.
    if (rhs.IO_ != 0) 
      IO_ = SetupFileIO( fileType_ );
    else
      IO_ = 0;
  }
  return *this;
}

// DESTRUCTOR
CpptrajFile::~CpptrajFile() {
   //fprintf(stderr,"CPPTRAJFILE DESTRUCTOR\n");
   CloseFile();
   if (IO_ != 0) delete IO_;
}
#ifdef MPI
/** Open the file using MPI file routines. */
int CpptrajFile::ParallelOpenFile(AccessType accessIn, Parallel::Comm const& commIn, bool sharedWrite)
{
  if (IO_ == 0) {
    mprinterr("Internal Error: CpptrajFile has not been set up.\n");
    return 1;
  }
  // This will currently only work for fileType_ STANDARD
  if (fileType_ != STANDARD) {
    mprinterr("Error: Parallel file access not supported for file type '%s'\n",
              FileTypeName[fileType_]);
    return 1;
  }
  // This will NOT work for streams.
  if (isStream_) {
    mprinterr("Error: Parallel file access not supported for streams.\n");
    return 1;
  }
  if (isOpen_) CloseFile();
  // TODO Save serial IO object?
  if (sharedWrite)
    fileType_ = MPISHARED;
  else
    fileType_ = MPIFILE;
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;
  ((FileIO_Mpi*)IO_)->SetComm( commIn );
  return OpenFile( accessIn );
}
#endif

// CpptrajFile::OpenFile()
/** Open the file. If already open, reopen.
  * NOTE: UPDATE access currently only used by Traj_CharmmDcd to update frame count.
  */
int CpptrajFile::OpenFile(AccessType accessIn) {
  if (IO_ == 0) {
    mprinterr("Internal Error: CpptrajFile has not been set up.\n");
    return 1;
  }
  int err = 0;
  if (isOpen_) CloseFile();
  if (isStream_) {
    switch (accessIn) {
      case READ : err = IO_->OpenStream( FileIO::STDIN ); break;
      case WRITE: err = IO_->OpenStream( FileIO::STDOUT); break;
      default:
        mprinterr("Internal Error: %s access not supported for file streams.\n",
                  AccessTypeName[accessIn]);
        err = 1;
    }
    if (debug_ > 0 && err == 0)
      rprintf("Opened stream %s\n", fname_.full());
  } else {
    if (fname_.empty()) {
      mprinterr("Internal Error: CpptrajFile file name is empty.\n");
      err = 1;
    } else {
      switch (accessIn) {
        case READ:   err = IO_->Open(fname_.full(), "rb"); break;
        case WRITE:  err = IO_->Open(fname_.full(), "wb"); break;
        case APPEND: err = IO_->Open(fname_.full(), "ab"); break;
        case UPDATE: err = IO_->Open(fname_.full(), "r+b"); break;
      }
      if (debug_ > 0 && err == 0)
        rprintf("Opened file %s with access %s\n", fname_.full(), AccessTypeName[accessIn]);
    }
  }
  if (err == 0)
    isOpen_ = true;
  else {
    if (debug_ > 0)
      rprinterr("Could not open %s with access %s\n", fname_.full(), AccessTypeName[accessIn]);
    mprinterr("Error: File '%s': %s\n", fname_.full(), strerror( errno ));
  }
  return err;
}

// CpptrajFile::CloseFile()
void CpptrajFile::CloseFile() {
  if (isOpen_) {
    IO_->Close();
    if (debug_>0) rprintf("Closed %s.\n", fname_.full());
    isOpen_=false;
#   ifdef MPI
    // Restore standard IO object.
    if (IsMPI()) {
      delete IO_;
      fileType_ = STANDARD;
      IO_ = SetupFileIO( fileType_ );
      if (IO_ == 0)
        mprinterr("Internal Error: Could not reset file '%s' from parallel to serial.\n",
                  fname_.full());
    }
#   endif
  }
}

// -----------------------------------------------------------------------------
// CpptrajFile::Printf()
/** Take the formatted string and write it to file using Write.
  */
void CpptrajFile::Printf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(linebuffer_,format,args);
  IO_->Write(linebuffer_, strlen(linebuffer_));
  va_end(args);
}

std::string CpptrajFile::GetLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE) != 0) {
    //mprinterr("Error: Getting line from %s\n", fname_.full());
    return std::string();
  }
  return std::string(linebuffer_);
}

const char* CpptrajFile::NextLine() {
  if (IO_->Gets(linebuffer_, BUF_SIZE) != 0) {
    //mprinterr("Error: Reading line from %s\n", fname_.full());
    return 0;
  }
  return linebuffer_;
}

// -----------------------------------------------------------------------------
// CpptrajFile::UncompressedSize()
off_t CpptrajFile::UncompressedSize() const {
  if (compressType_ == NO_COMPRESSION)
    return file_size_;
  else
    return uncompressed_size_;
}

// CpptrajFile::Reset()
/** Close file if open, reset all file information.
  */
void CpptrajFile::Reset() {
  CloseFile();
  if (IO_!=0) delete IO_;
  IO_ = 0;
  fname_.clear();
  isOpen_ = false;
  isStream_ = false;
  uncompressed_size_ = 0UL;
  compressType_ = NO_COMPRESSION;
  isDos_ = 0;
}

// -----------------------------------------------------------------------------
// CpptrajFile::OpenRead()
int CpptrajFile::OpenRead(FileName const& nameIn) {
  if (SetupRead( nameIn, debug_ )) return 1; 
  return OpenFile();
}

// CpptrajFile::OpenWrite()
int CpptrajFile::OpenWrite(FileName const& nameIn) {
  if (SetupWrite(nameIn, debug_)) return 1;
  return OpenFile();
}

// CpptrajFile::OpenWriteNumbered()
// NOTE: File MUST be previously set up. Primarily for use with traj files.
int CpptrajFile::OpenWriteNumbered(int numIn, bool prepend) {
  if (isStream_) {
    mprinterr("Internal Error: CpptrajFile::OpenWriteNumbered cannot be used with streams.\n");
    return 1;
  }
  if (prepend) {
    FileName newName = fname_.PrependExt( "." + integerToString(numIn) );
    if (IO_->Open( newName.full(), "wb")) return 1;
  } else {
    std::string newName = AppendNumber( fname_.Full(), numIn );
    if (IO_->Open( newName.c_str(), "wb")) return 1;
  }
  isOpen_ = true;
  return 0;
}

// CpptrajFile::OpenAppend()
int CpptrajFile::OpenAppend(FileName const& nameIn) {
  if (SetupAppend(nameIn, debug_)) return 1;
  return OpenFile();
}

// -----------------------------------------------------------------------------
// CpptrajFile::SetupRead()
/** Set up file for reading. Will autodetect the type.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupRead(FileName const& nameIn, int debugIn) {
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = READ;
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for READ.\n", nameIn.full());
  // If nameIn is empty assume reading from STDIN desired. 
  if (nameIn.empty()) {
    isStream_ = true;
    // file type must be STANDARD for streams
    fileType_ = STANDARD;
    fname_.SetFileName_NoExpansion("STDIN");
    IO_ = SetupFileIO( fileType_ );
  } else {
    isStream_ = false;
    // Check if file exists. If not, fail silently
    if (!File::Exists( nameIn )) return 1;
    fileType_ = UNKNOWN_TYPE;
    // Determine file type. This sets up IO and determines compression. 
    if (ID_Type( nameIn.full() )) return 1;
    // Set up filename; sets base filename and extensions
    fname_ = nameIn;
  } 
  if (debug_>0)
    rprintf("\t[%s] is type %s with access READ\n", fname_.full(), FileTypeName[fileType_]);
  return 0;
}

// CpptrajFile::SetupWrite()
int CpptrajFile::SetupWrite(FileName const& nameIn, int debugIn) {
  return SetupWrite(nameIn, UNKNOWN_TYPE, debugIn);
}

// CpptrajFile::SetupWrite()
/** Set up file for writing with the given type. If no filename is given 
  * this indicates STDOUT. If no type is specified attempt to detect
  * from the compression extension.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupWrite(FileName const& nameIn, FileType typeIn, int debugIn) 
{
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = WRITE;
  fileType_ = typeIn;
  // If nameIn is empty assume writing to STDOUT desired.
  if (nameIn.empty()) {
    isStream_ = true;
    // file type must be STANDARD for streams
    fileType_ = STANDARD;
    fname_.SetFileName_NoExpansion("STDOUT");
  } else {
    isStream_ = false;
    // Set up filename; sets base filename and extension
    fname_ = nameIn;
  }
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for WRITE.\n", fname_.full());
  // If file type is not specified, try to determine from filename extension
  if (fileType_ == UNKNOWN_TYPE) {
    if (fname_.Compress() == ".gz")
      fileType_ = GZIPFILE;
    else if (fname_.Compress() == ".bz2")
      fileType_ = BZIP2FILE;
    else
      fileType_ = STANDARD;
  }
  // Setup IO based on type.
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;
  if (debug_>0)
    rprintf("\t[%s] is type %s with access WRITE\n", fname_.full(), FileTypeName[fileType_]);
  return 0;
}

// CpptrajFile::SetupAppend()
/** Set up the file for appending. Will first set up for read to determine
  * the type and format. Set up for write if file does not exist.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupAppend(FileName const& nameIn, int debugIn) {
  // Make append to null an error
  if (nameIn.empty()) {
    mprinterr("Error: SetupAppend(): No filename specified\n");
    return 1;
  }
  // NOTE: File will be cleared and debug set by either SetupRead/SetupWrite
  if (File::Exists(nameIn)) {
    // If file exists, first set up for read to determine type and format.
    if (SetupRead(nameIn, debugIn)!=0) return 1;
    access_ = APPEND;
  } else {
    // File does not exist, just set up for write.
    if (SetupWrite(nameIn, debugIn)!=0) return 1;
    if (debug_>0)
      mprintf("Warning: %s not accessible, changed access from APPEND to WRITE.\n",
              fname_.full());
  }
  // Appending and compression not supported.
  if (IsCompressed()) {
    mprinterr("Error: Appending to compressed files is not supported.\n");
    return 1;
  }
  if (debug_>0)
    rprintf("\t[%s] is type %s with access APPEND\n", fname_.full(), FileTypeName[fileType_]);
  return 0;
}

// -----------------------------------------------------------------------------
// CpptrajFile::SetupFileIO()
/** Set up the IO based on given file type. */
FileIO* CpptrajFile::SetupFileIO(FileType typeIn) {
  switch (typeIn) {
    case STANDARD  : return (new FileIO_Std());
    case GZIPFILE  : 
#ifdef HASGZ
      return new FileIO_Gzip(); 
#else
      mprinterr("Error: Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 0;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      return (new FileIO_Bzip2());
#else
      mprinterr("Error: Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 0;
#endif
    break;
#ifdef MPI
    case MPIFILE   : return (new FileIO_Mpi());
    case MPISHARED : return (new FileIO_MpiShared());
#else
    case MPIFILE   :
    case MPISHARED :
      mprinterr("Error: Compiled without MPI support. Recompile with -DMPI\n");
      return 0;
#endif
      break;
    //case ZIPFILE   : return (new ZipFile()); break;
    default : 
      mprinterr("Error: Unrecognized file type.\n");
  }
  return 0;
}

// CpptrajFile::ID_Type() 
/** Attempt to identify the file type for filenameIn. Also set file_size,
  * uncompressed_size, and compressType.
  * FIXME: Will have to be modified to use OpenFile if STDIN ever enabled.
  */
int CpptrajFile::ID_Type(const char* filenameIn) {
  if (filenameIn == 0) return 1;
  // Get basic file information
  struct stat frame_stat;
  if (stat(filenameIn, &frame_stat) == -1) {
    mprinterr( "Error: Could not find file status for %s\n", filenameIn);
    if (debug_>0) 
      perror("     Error from stat: ");
    return 1;
  }
  file_size_ = frame_stat.st_size;
  // Start off every file as a standard file
  fileType_ = STANDARD;
  IO_ = new FileIO_Std();
  // ID by magic number - open for binary read access
  if ( IO_->Open(filenameIn, "rb") ) { 
    mprintf("Could not open %s for hex signature read.\n", filenameIn);
    return 1;
  }
  // Read first 3 bytes
  unsigned char magic[3];
  magic[0] = 0; 
  magic[1] = 0; 
  magic[2] = 0;
  IO_->Read(magic, 3);
  IO_->Close();
  if (debug_>0) mprintf("\t    Hex sig: %x %x %x", magic[0],magic[1],magic[2]);
  // Check compression
  if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8)) {
    if (debug_>0) mprintf(", Gzip file.\n");
    compressType_ = GZIP;
    fileType_ = GZIPFILE;
  } else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68)) {
    if (debug_>0) mprintf(", Bzip2 file.\n");
    compressType_ = BZIP2;
    fileType_ = BZIP2FILE;
  } else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3)) {
    if (debug_>0) mprintf(", Zip file.\n");
    compressType_ = ZIP;
    fileType_ = ZIPFILE;
  } else {
    if (debug_>0) mprintf(", No compression.\n");
  }
  // Assign the appropriate IO type based on the file type
  delete (FileIO_Std*) IO_;
  IO_ = SetupFileIO( fileType_ );
  if (IO_ == 0) return 1;

  // If the file is compressed, get the uncompressed size
  // For standard files this just returns 0UL
  // Standard file size is in the frame_stat struct
  uncompressed_size_ = IO_->Size(filenameIn);

  // Additional file characteristics
  linebuffer_[0]='\0';
  if (IO_->Open(filenameIn, "rb")!=0) return 1; 
  IO_->Gets(linebuffer_,BUF_SIZE);
  IO_->Close();

  // Check for terminal CR before newline, indicates DOS file
  size_t i = strlen(linebuffer_);
  if ( i>1 ) {
    if (linebuffer_[ i - 2 ] == '\r') {
      if (debug_>0) mprintf("  [DOS]");
      isDos_ = 1;
    }
  }
  return 0;
}
