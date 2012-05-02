// CpptrajFile
// TODO: Replace sprintf/atof with sstream functs?
#include <cstring> // strlen 
#include <sys/stat.h> // stat
#include <cstdio> // sprintf, vsprintf
#include <cstdlib> // atof
#include <cstdarg> // va_X functions
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "MpiRoutines.h" // For worldrank in Rank_printf
#  include "FileIO_Mpi.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

const char CpptrajFile::FileTypeName[6][13] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};

const char CpptrajFile::AccessTypeName[3][2] = {
  "R", "W", "A"
};

const size_t CpptrajFile::BUF_SIZE = 83;

// CONSTRUCTOR
CpptrajFile::CpptrajFile() :
  IO(NULL),
  access_(READ),
  isDos_(0),
  uncompressed_size_(0UL),
  file_size_(0UL),
  compressType_(NO_COMPRESSION),
  debug_(0),
  isOpen_(false),
  fileType_(STANDARD),
  filename_(NULL)
{}

// Copy Constructor
CpptrajFile::CpptrajFile(const CpptrajFile &rhs) {
  // Even if file is open, copy it closed?
  isOpen_ = false;
  debug_ = rhs.debug_;
  access_ = rhs.access_;
  uncompressed_size_ = rhs.uncompressed_size_;
  file_size_ = rhs.file_size_;
  fileType_ = rhs.fileType_;
  FileName_ = rhs.FileName_;
  if (FileName_.empty())
    filename_=NULL;
  else
    filename_ = (char*)FileName_.c_str();
  basefilename_ = rhs.basefilename_;
  Ext_ = rhs.Ext_;
  compressType_ = rhs.compressType_;
  isDos_ = rhs.isDos_;
  // If charbuffer, assign here
  // Set up the IO object
  // NOTE: Should probably throw an exception if this fails.
  if (rhs.IO!=NULL) 
    SetupFileIO();
  else
    IO = NULL;
}

// Assignment
CpptrajFile &CpptrajFile::operator=(const CpptrajFile &rhs) {
  // Self-assignment
  if (this == &rhs) return *this;
  // Deallocate
  CloseFile();
  if (IO!=NULL) delete IO;
  // Allocate and copy
  debug_ = rhs.debug_;
  access_ = rhs.access_;
  uncompressed_size_ = rhs.uncompressed_size_;
  file_size_ = rhs.file_size_;
  fileType_ = rhs.fileType_;
  FileName_ = rhs.FileName_;
  if (FileName_.empty())
    filename_=NULL;
  else
    filename_ = (char*)FileName_.c_str();
  basefilename_ = rhs.basefilename_;
  Ext_ = rhs.Ext_;
  compressType_ = rhs.compressType_;
  isDos_ = rhs.isDos_;
  // If charbuffer, assign here
  // Set up the IO object
  // NOTE: Should probably throw an exception if this fails.
  if (rhs.IO!=NULL) 
    SetupFileIO();
  else
    IO = NULL;
  return *this;
}

// DESTRUCTOR
CpptrajFile::~CpptrajFile() {
   //fprintf(stderr,"CPPTRAJFILE DESTRUCTOR\n");
   CloseFile();
   if (IO!=NULL) delete IO;
}

// CpptrajFile::IsOpen()
bool CpptrajFile::IsOpen() { 
  if (isOpen_) return true; 
  return false; 
}

// CpptrajFile::Name()
const char *CpptrajFile::Name() {
  return FileName_.c_str();
}

// CpptrajFile::FullPathName()
std::string CpptrajFile::FullPathName() {
  return FileName_;
}

// CpptrajFile::BaseName()
const char *CpptrajFile::BaseName() {
  return basefilename_.c_str();
}

// CpptrajFile::Extension()
std::string CpptrajFile::Extension() {
  return Ext_;
}

// CpptrajFile::IsDos()
bool CpptrajFile::IsDos() {
  if (isDos_==1) return true;
  return false;
}

// CpptrajFile::IsCompressed()
bool CpptrajFile::IsCompressed() {
  if (compressType_ == NO_COMPRESSION) return false;
  return true;
}

// CpptrajFile::OpenFile()
/** Open the file. If already open, reopen.  */
int CpptrajFile::OpenFile() {
  if (isOpen_) CloseFile();

  switch (access_) {
    case READ:
      if (filename_ == NULL) {
        mprinterr("Error: CpptrajFile: Filename is NULL.\n");
        return 1;
      }
      if ( IO->Open(filename_, "rb")  ) { 
        rprintf("Could not open %s for reading.\n",filename_);
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for reading.\n",filename_);
      break;
    case APPEND:
      if (filename_ == NULL) {
        mprinterr("Error: CpptrajFile: Filename is NULL.\n");
        return 1;
      }
      if ( IO->Open(filename_, "ab") ) {
        rprintf("Could not open %s for appending.\n",filename_);
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for appending.\n",filename_);
      break;
    case WRITE:
      if ( IO->Open(filename_, "wb") ) { 
        rprintf("Could not open %s for writing.\n",filename_);
        return 1;
      }
      if (debug_>0) rprintf("Opened %s for writing.\n",filename_);
      break;
  }
      
  isOpen_=true;
  return 0;
}

// CpptrajFile::CloseFile()
/** Close the file.  */
void CpptrajFile::CloseFile() {
  if (isOpen_) {
    IO->Close();
    if (debug_>0) rprintf("Closed %s.\n",filename_);
    isOpen_=false;
  }
}

// -----------------------------------------------------------------------------
// CpptrajFile::Printf()
/** Take the formatted string and write it to file using Write.
  */
void CpptrajFile::Printf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(printf_buffer_,format,args);
  IO->Write(printf_buffer_, strlen(printf_buffer_));
  va_end(args);
}

// CpptrajFile::Rank_printf()
/** When MPI, printf only for the specified rank. If no MPI, behaves just
  * like above Printf.
  */
void CpptrajFile::Rank_printf(int rank, const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsprintf(printf_buffer_,format,args);
#ifdef MPI
    if (worldrank==rank)
#endif
      IO->Write(printf_buffer_, strlen(printf_buffer_));
  va_end(args);
}

// -----------------------------------------------------------------------------
// CpptrajFile::SetBaseFilename()
/** Set filename with full path. Strip leading path from input filename to 
  * determine the base filename. Also determine the file extension.
  */
void CpptrajFile::SetBaseFilename(const char *filenameIn) {
  //mprintf("DEBUG: Called SetBaseFilename with [%s]\n",filenameIn);
  // NULL filename allowed for WRITE (indicates STDOUT)
  if (filenameIn==NULL) {
    filename_ = NULL;
    return;
  }
  // Assign filename with full path
  FileName_.assign(filenameIn);
  filename_ = (char*)FileName_.c_str();

  // Get position of last occurence of '/' to determine base filename
  size_t found = FileName_.find_last_of("/");
  if (found == std::string::npos)
    basefilename_ = FileName_;
  else
    basefilename_ = FileName_.substr(found+1);

  if (debug_>0) {
    mprintf("\tSetBaseFilename: Filename is %s\n",filename_);
    mprintf("\t                 Base filename is %s\n",basefilename_.c_str());
  }

  // Get the filename extension
  found = basefilename_.find_last_of(".");
  if (found == std::string::npos) {
    Ext_.clear();
    if (debug_>0) mprintf("\t                 No extension.\n");
  } else {
    Ext_ = basefilename_.substr(found);
    if (debug_>0)
      mprintf("\t                 Ext= %s  Len= %lu\n",Ext_.c_str(),Ext_.size());
  }
}

// CpptrajFile::Reset()
/** Close file if open, reset all file information.
  */
void CpptrajFile::Reset() {
  CloseFile();
  if (IO!=NULL) delete IO;
  IO = NULL;
  FileName_.clear();
  filename_ = NULL;
  basefilename_.clear();
  Ext_.clear();
  isOpen_ = false;
  uncompressed_size_ = 0UL;
  compressType_ = NO_COMPRESSION;
  isDos_ = 0;
}

// CpptrajFile::SetupRead()
/** Set up file for reading. Will autodetect the type and format.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupRead(const char *filenameIn, int debugIn) {
  // NULL filename not allowed
  if (filenameIn==NULL) {
    mprinterr("Internal Error: NULL filename specified for READ.\n");
    return 1;
  }
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = READ;
  // Setting fileFormat to unknown implicitly requests detection.
  //fileFormat_ = UNKNOWN_FORMAT;
  fileType_ = UNKNOWN_TYPE;
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for READ.\n",filenameIn);
  // Perform tilde-expansion
  std::string expandedName = tildeExpansion( filenameIn );
  if (expandedName.empty()) {
    // Failure of tilde-expansion probably means file does not exist.
    //mprinterr("Interal Error: CpptrajFile: Tilde-expansion failed.\n");
    return 1;
  }
  // Set up filename; sets base filename and extension
  // TODO: Pass in string
  SetBaseFilename((char*)expandedName.c_str());
  // Determine file type. This sets up IO. 
  if (ID_Type()) return 1;
  if (debug_>0)
    rprintf("\t[%s] is type %s with access READ\n",filename_,FileTypeName[fileType_]);
    //rprintf("\t[%s] is format %s and type %s with access READ\n",filename_,
    //        File_Format(fileFormat_), File_Type(fileType_));
  return 0;
}

// CpptrajFile::SetupAppend()
/** Set up the file for appending. Will first set up for read to determine
  * the type and format.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupAppend(const char *filenameIn, int debugIn) {
  // First set up for read to determine type and format.
  if (SetupRead(filenameIn,debugIn)!=0) return 1;
  access_ = APPEND;
  if (debug_>0)
    mprintf("CpptrajFile: Changed %s access to APPEND.\n",filename_);
  return 0;
}

// SetupWrite()
/** Set up file for writing with the given format and type. If a NULL filename
  * is given this indicates STDOUT.
  * \return 0 on success, 1 on error.
  */
//int CpptrajFile::SetupWrite(char *filenameIn, FileFormat fmtIn, FileType typeIn, int debugIn)
int CpptrajFile::SetupWrite(const char *filenameIn, int debugIn) 
{
  // Clear file, set debug level
  Reset();
  debug_ = debugIn;
  access_ = WRITE;
  if (debug_>0)
    mprintf("CpptrajFile: Setting up %s for WRITE.\n",filenameIn);
  // Set up filename; sets base filename and extension
  SetBaseFilename(filenameIn);
  // If file type is not specified, try to determine from filename extension
  //fileType_ = DetermineType(UNKNOWN_TYPE, (char*)Ext_.c_str());
  if (Ext_ == ".gz")
    fileType_ = GZIPFILE;
  else if (Ext_ == ".bz2")
    fileType_ = BZIP2FILE;
  else
    fileType_ = STANDARD; 
  // If file format not specified, try to determine from filename extension
  //fileFormat_ = DetermineFormat(fmtIn, (char*)Ext_.c_str());
  // Setup IO based on type.
  if (SetupFileIO()) return 1;
  if (debug_>0)
    rprintf("\t[%s] is type %s with access WRITE\n",filename_,FileTypeName[fileType_]);
    //rprintf("\t[%s] is format %s and type %s with access WRITE\n",filename_,
    //        File_Format(fileFormat_), File_Type(fileType_));
  return 0;
}

// Cpptraj::SetupFileIO()
/** Set up the IO based on previously determined file type.
  */
int CpptrajFile::SetupFileIO() {
  switch (fileType_) {
    case STANDARD  : IO = new FileIO_Std();  break;
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new FileIO_Gzip(); 
#else
      mprinterr("Error: (%s): Compiled without Gzip support. Recompile with -DHASGZ\n",
                filename_);
      return 1;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      IO = new FileIO_Bzip2();
#else
      mprinterr("Error: (%s): Compiled without Bzip2 support. Recompile with -DHASBZ2\n",
                filename_);
      return 1;
#endif
    break;
    case MPIFILE   : 
#ifdef MPI
      IO = new FileIO_Mpi();
#else
      mprinterr("Error: (%s): Compiled without MPI support. Recompile with -DMPI\n",
                filename_);
      return 1;
#endif
      break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    default : 
      mprinterr("Error: (%s): Unrecognized file type.\n",filename_);
      return 1;
  }
  if (IO==NULL) return 1;
  return 0;
}

// CpptrajFile::ID_Type() 
/** Open the file specified by filename for READ or APPEND access. Attempt to 
  * identify the file type.
  */
int CpptrajFile::ID_Type() {
  unsigned char magic[3];
  char buffer1[BUF_SIZE];
  struct stat frame_stat;

  //mprintf("DEBUG: ID_Type( %s )\n",filename_);
  // Get basic file information
  // An error here means file probably doesnt exist. Dont print an error at 
  // basic debug level since this could also be used to test if file exists.
  if (stat(filename_, &frame_stat) == -1) {
    mprinterr( "Error: Could not find file status for %s\n", filename_);
    if (debug_>0) {
      perror("     Error from stat: ");
    }
    return 1;
  }
  file_size_ = frame_stat.st_size;

  // Start off every file as a standard file - may need to change for MPI
  fileType_ = STANDARD;
  IO = new FileIO_Std();

  // ID by magic number - open for binary read access
  if ( IO->Open(filename_, "rb") ) { 
    mprintf("Could not open %s for hex signature read.\n",filename_);
    return 1;
  }

  // Read first 3 bytes
  magic[0] = 0; magic[1] = 0; magic[2] = 0;
  IO->Read(magic  ,1);
  IO->Read(magic+1,1);
  IO->Read(magic+2,1);
  IO->Close();
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

  // Appending and compression not supported.
  if (access_==APPEND && compressType_!=NO_COMPRESSION) {
    mprintf("Error: Appending to compressed files is not supported.\n");
    return 1;
  }

  // Assign the appropriate IO type based on the file type
  delete (FileIO_Std*) IO;
  IO = NULL;
  if (SetupFileIO()) return 1;

  // If the file is compressed, get the uncompressed size
  // For standard files this just returns 0UL
  // Standard file size is in the frame_stat struct
  uncompressed_size_ = IO->Size(filename_);

  // Additional file characteristics
  buffer1[0]='\0';
  if (IO->Open(filename_,"rb")!=0) return 1; 
  IO->Gets(buffer1,BUF_SIZE);
  IO->Close();

  // Check for terminal CR before newline, indicates DOS file
  size_t i = strlen(buffer1);
  if ( i>1 ) {
    if (buffer1[ i - 2 ] == '\r') {
      if (debug_>0) mprintf("  [DOS]");
      isDos_ = 1;
    }
  }
  return 0;
}

