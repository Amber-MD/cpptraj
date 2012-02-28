// CpptrajFile
/* Compiler Defines:
 * - USE_CHARBUFFER: Use CharBuffer to buffer entire file
 */ 
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <sys/stat.h>
#include "CpptrajFile.h"
#include "NetcdfRoutines.h"
#include "PDBfileRoutines.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// File Types
#include "FileIO_Std.h"
#ifdef HASGZ
#  include "FileIO_Gzip.h"
#endif
#ifdef MPI
#  include "FileIO_Mpi.h"
#endif
#ifdef HASBZ2
#  include "FileIO_Bzip2.h"
#endif

// CONSTRUCTOR
CpptrajFile::CpptrajFile() {
  IO=NULL;
  isOpen=false;
  uncompressed_size=0UL;
  compressType=NO_COMPRESSION;
  access=READ;
  debug=0;
  fileFormat=UNKNOWN_FORMAT;
  fileType=STANDARD;
  filename=NULL;
  basefilename=NULL;
  Ext=NULL;
  isDos=0;
}

// DESTRUCTOR
CpptrajFile::~CpptrajFile() {
   //fprintf(stderr,"CPPTRAJFILE DESTRUCTOR\n");
   CloseFile();
   if (IO!=NULL) delete IO;
   if (filename!=NULL) free(filename);
   if (basefilename!=NULL) free(basefilename);
   if (Ext!=NULL) free(Ext);
}

// CpptrajFile::OpenFile()
/** Open the file. If already open, reopen.  */
int CpptrajFile::OpenFile() {
  if (isOpen) CloseFile();

  switch (access) {
    case READ:
      if ( IO->Open(filename, "rb")  ) { // NOTE: use rb as mode instead?
        rprintf("Could not open %s for reading.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for reading.\n",filename);
      break;
    case APPEND:
      if ( IO->Open(filename, "ab") ) {
        rprintf("Could not open %s for appending.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for appending.\n",filename);
      break;
    case WRITE:
      if ( IO->Open(filename, "wb") ) { // NOTE: Use wb as mode?
        rprintf("Could not open %s for writing.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for writing.\n",filename);
      break;
  }
      
  isOpen=true;
  return 0;
}

// CpptrajFile::CloseFile()
/** Close the file.  */
void CpptrajFile::CloseFile() {
  if (isOpen) {
    IO->Close();
    if (debug>0) rprintf("Closed %s.\n",filename);
    isOpen=false;
  }
}

#ifdef USE_CHARBUFFER
// ---------- BUFFERED FILE ROUTINES -------------------------------------------
// CpptrajFile::OpenFileBuffered()
/** Open the file. Buffer with CharBuffer. */
int CpptrajFile::OpenFileBuffered() {
  size_t blockSize = BUFFER_SIZE;
  if (isOpen) CloseFile();

  switch (access) {
    case READ:
      if ( IO->Open(filename, "rb")  ) { // NOTE: use rb as mode instead?
        rprintf("Could not open %s for reading.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for reading.\n",filename);
      // Allocate char buffer to read in entire file if possible.
      if (compressType == NO_COMPRESSION || uncompressed_size == 0)
        blockSize = file_size;
      else
        blockSize = uncompressed_size;
      c_buffer.Allocate( blockSize );
      // Read in entire file
      if (IO->Read(c_buffer.BufferPtr(), sizeof(char), blockSize)==-1) return 1;
      // If compressed file and uncompressed size not known keep reading in
      // blocks until EOF.
      if (compressType!=NO_COMPRESSION && uncompressed_size == 0) {
        bool readMore = true;
        while (readMore) {
          c_buffer.IncreaseSize( blockSize );
          int numread = IO->Read(c_buffer.BufferPtr(), sizeof(char), blockSize);
          if (numread==-1 || numread < (int)blockSize) readMore=false;
          // NOTE: Resize buffer here to match actual # chars read in?
        }
      }
      break;
    case APPEND:
/*      if ( IO->Open(filename, "ab") ) {
        rprintf("Could not open %s for appending.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for appending.\n",filename);
      // Allocate char buffer with default size
      c_buffer.Allocate( BUFFER_SIZE );*/
      mprinterr("Internal Error: Buffered writes currently disabled.\n");
      return 1;
      break;
    case WRITE:
/*      if ( IO->Open(filename, "wb") ) { // NOTE: Use wb as mode?
        rprintf("Could not open %s for writing.\n",filename);
        return 1;
      }
      if (debug>0) rprintf("Opened %s for writing.\n",filename);
      // Allocate char buffer with default size
      c_buffer.Allocate( BUFFER_SIZE );*/
      mprinterr("Internal Error: Buffered writes currently disabled.\n");
      return 1;
      break;
  }

  isOpen=true;
  return 0;
}

// CpptrajFile::ReadBuffered()
/** Read BufferSize bytes into the CharBuffer. */
/*int CpptrajFile::ReadBuffered() {
  char *bufferptr = c_buffer.BufferPtr();
  if (IO->Read(bufferptr, sizeof(char), c_buffer.BufferSize())==-1) return 1;
  return 0;
}*/

// CpptrajFile::Gets()
/// Gets from the CharBuffer.
int CpptrajFile::Gets(char *str, int num) {
  return c_buffer.Gets(str,num);
}

// CpptrajFile::Rewind()
/// Rewind the CharBuffer.
void CpptrajFile::Rewind() {
  c_buffer.Rewind();
}

// CpptrajFile::Read()
/// Read from the CharBuffer.
int CpptrajFile::Read(void *str, size_t numbytes) {
  return c_buffer.Read(str,numbytes);
}
// -----------------------------------------------------------------------------
#endif

// CpptrajFile::SetBaseFilename()
/** Strip leading path from input filename. Use strtok routine to separate 
  * filename by / and use the last string as the base filename. Internal 
  * filename is not used since strtok modifies the char array.
  * Also determine the file extension.
  */
void CpptrajFile::SetBaseFilename() {
  char *ptr, *tempFilename;
  int i, nameLen;

  if (filename==NULL) return;
  tempFilename=(char*) malloc( (strlen(filename)+1) * sizeof(char));
  strcpy(tempFilename, filename);

  basefilename=NULL;
  ptr = strtok(tempFilename,"/");
  while (ptr != NULL) {
    //printf ("DEBUG: %s\n",ptr);
    basefilename = (char*) realloc(basefilename, (strlen(ptr)+1) * sizeof(char));
    strcpy(basefilename,ptr);
    ptr = strtok(NULL,"/");
  }
  if (debug>0) {
    mprintf("\tSetBaseFilename: Filename is %s\n",filename);
    mprintf("\t                 Base filename is %s\n",basefilename);
  }
  free(tempFilename);

  // Get the filename extension
  nameLen=strlen(basefilename);
  for (i=nameLen-1; i>=0; i--)
    if (basefilename[i]=='.') break;
  // No extension
  if (i<0) {
    Ext=NULL;
    if (debug>0) 
      mprintf("\t                 No extension.\n");
    return;
  } 
  Ext = (char*) malloc( (strlen(basefilename+i) + 1) * sizeof(char));
  strcpy(Ext, basefilename+i);
  if (debug>0)
    mprintf("\t                 Ext= %s  Len= %lu\n",Ext,strlen(basefilename+i));
}

// CpptrajFile::SetupFile()
/** Set up the given file for the specified access, assume DATAFILE format
  * (no autodetection of format). Implicitly requests autodetection of the 
  * file type.
  */
int CpptrajFile::SetupFile(char *filenameIn, AccessType accessIn, int debugIn) {
  return SetupFile(filenameIn,accessIn,DATAFILE,UNKNOWN_TYPE,debugIn);
}

// CpptrajFile::SetupFile()
/** Sets the file name, access type (R/W/A), file type and file format (for 
  * WRITE), and debug level. If called with READ or append file type and format
  * will be determined by SetupRead. If called with WRITE the given type and 
  * format will be used; if called with UNKNOWN type and format they will
  * try to be determined by the given file extension. 
  * Can be called with NULL for write, this will write to stdout.
  * \return 0 on success, 1 on error.
  */
int CpptrajFile::SetupFile(char *filenameIn, AccessType accessIn, 
                         FileFormat fileFormatIn, FileType fileTypeIn, 
                         int debugIn) 
{
  // If file has been previously setup / opened, clear file.
  CloseFile();
  if (IO!=NULL) delete IO;
  if (filename!=NULL) free(filename);
  if (basefilename!=NULL) free(basefilename);
  if (Ext!=NULL) free(Ext);
  IO=NULL;
  filename=NULL;
  basefilename=NULL;
  Ext=NULL;
  isOpen=false;
  uncompressed_size=0UL;
  compressType=NO_COMPRESSION;
  isDos=0;
  debug=debugIn;
  if (debug>0) {
    mprintf("CpptrajFile: [%s] FMT %s, TYPE %s, ACC %s\n",filenameIn,
           File_Format(fileFormatIn),File_Type(fileTypeIn),
           File_Access(accessIn));
  }
  // Store tilde-expanded filename on read/append
  if (accessIn!=WRITE) {
    filename = tildeExpansion(filenameIn,debug);
    if (filename==NULL) {
      mprinterr("Error: CpptrajFile::SetupFile: NULL filename specified on READ or APPEND\n");
      return 1;
    }
  // On WRITE, just store filename if not NULL (NULL indicates STDOUT)
  } else { 
    if (filenameIn!=NULL) {
      filename=(char*) malloc( (strlen(filenameIn)+1) * sizeof(char) );
      strcpy(filename,filenameIn);
    }
  }
  // Store base filename and determine filename extension
  this->SetBaseFilename();

  access=accessIn;
  switch (access) {
    case APPEND:
      // If appending, type and format MUST be determined.
      if ( SetupRead()  ) return 1;
      break;
    case READ:  
      // SetupRead will determine format if fileFormatIn is UNKNOWN.
      // SetupRead always determines type.
      fileFormat = fileFormatIn;
      if ( SetupRead()  ) return 1; 
      break;
    case WRITE: 
      // If file type is not specified, try to determine from filename
      // NOTE: Also determine format from name??
      fileType=DetermineType(fileTypeIn,Ext);
      fileFormat=DetermineFormat(fileFormatIn,Ext);
      if ( SetupWrite() ) return 1; 
      break;
    default: return 1;
  }
  if (debug>0)
    rprintf("\t[%s] is format %s and type %s with access %s\n",filename,
           File_Format(fileFormat),File_Type(fileType),File_Access(access));
  return 0;
}

// CpptrajFile::SetupWrite()
/** Set up file with specified type for writing. fileFormat is set by SetupFile.
  * NOTE: Determine if compression is requested, either by arg or from name
  */
int CpptrajFile::SetupWrite() {

  if (debug>1) mprintf("DEBUG: Setting up WRITE for file %s\n",filename);
  // Eventually allow other file types
  //fileType=STANDARD;
  switch (fileType) {
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new FileIO_Gzip(); 
#else
      mprintf("Error: SetupWrite(%s):\n",filename);
      mprintf("       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      IO = new FileIO_Bzip2();
#else
      mprintf("Error: SetupWrite(%s):\n",filename);
      mprintf("       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#endif
    break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    case STANDARD  : IO = new FileIO_Std();  break;
    case MPIFILE   : 
#ifdef MPI
      IO = new FileIO_Mpi();
#else
      mprintf("Error: SetupWrite(%s):\n",filename);
      mprintf("       Compiled without MPI support. Recompile with -DMPI\n");
      return 1;
#endif
      break;
    default : 
      mprintf("CpptrajFile::SetupWrite: Unrecognized file type.\n");
      return 1;
  }

  return 0;
}

// CpptrajFile::SetupRead() 
/** Open the file specified by filename for READ or APPEND access. Attempt to 
  * identify the file type and format.
  */
int CpptrajFile::SetupRead() {
  unsigned char magic[3];
  char buffer1[BUFFER_SIZE], buffer2[BUFFER_SIZE];
  char *CheckConventions; // Only used to check if netcdf is traj or restart
  float TrajCoord[10];
  int iamber[12];
  int i;
  struct stat frame_stat;

  //mprintf("DEBUG: Setting up read file %s\n",filename);
  // Get basic file information
  // An error here means file probably doesnt exist. Dont print an error at 
  // basic debug level since this could also be used to test if file exists.
  if (stat(filename, &frame_stat) == -1) {
    mprinterr( "Error: CpptrajFile::SetupRead: Could not find file status for %s\n", filename);
    if (debug>0) {
      perror("     Error from stat: ");
    }
    return 1;
  }
  file_size = frame_stat.st_size;

  // Start off every file as a standard file - may need to change for MPI
  IO = new FileIO_Std();

  // ID by magic number - open for binary read access
  if ( IO->Open(filename, "rb") ) { 
    mprintf("Could not open %s for hex signature read.\n",filename);
    return 1;
  }

  // Read first 3 bytes
  memset(magic,0,3*sizeof(unsigned char));
  IO->Read(magic  ,1,1);
  IO->Read(magic+1,1,1);
  IO->Read(magic+2,1,1);
  IO->Close();
  if (debug>0) mprintf("\t    Hex sig: %x %x %x", magic[0],magic[1],magic[2]);

  // Check compression
  if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8)) {
    if (debug>0) mprintf(", Gzip file.\n");
    compressType=GZIP;
    fileType=GZIPFILE;
  } else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68)) {
    if (debug>0) mprintf(", Bzip2 file.\n");
    compressType=BZIP2;
    fileType=BZIP2FILE;
  } else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3)) {
    if (debug>0) mprintf(", Zip file.\n");
    compressType=ZIP;
    fileType=ZIPFILE;
  } else {
    if (debug>0) mprintf(", No compression.\n");
  }

  // Appending and compression not supported.
  if (access==APPEND && compressType!=NO_COMPRESSION) {
    mprintf("Error: Appending to compressed files is not supported.\n");
    return 1;
  }

  // Assign the appropriate IO type based on the file
  delete (FileIO_Std*) IO;
  IO = NULL;
  switch (fileType) {
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new FileIO_Gzip(); 
#else
      mprintf("Error: SetupRead(%s):\n",filename);
      mprintf("       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE : 
#ifdef HASBZ2
      IO = new FileIO_Bzip2(); 
#else
      mprintf("Error: SetupRead(%s):\n",filename);
      mprintf("       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#endif
      break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    default        : IO = new FileIO_Std();  break;
  }

  // If the file is compressed, get the uncompressed size
  // For standard files this just returns 0UL
  // Standard file size is in the frame_stat struct
  uncompressed_size = IO->Size(filename);

  // If file format already determined, exit.
  if (fileFormat!=UNKNOWN_FORMAT) {
    if (debug>0) mprintf("CpptrajFile: %s: Not detecting file format.\n",filename);
    return 0;
  }

  // ========== Determine format if UNKNOWN ==========
  // NOTE: Should this section be moved to FileRoutines?
  // ---------- Read first 3 bytes again to determine format by magic number ----------
  IO->Open(filename,"rb"); // NOTE: Err Check
  memset(magic,0,3*sizeof(unsigned char));
  IO->Read(magic  ,1,1);
  IO->Read(magic+1,1,1);
  IO->Read(magic+2,1,1);
  IO->Close();
  if (debug>0) mprintf("\t    Hex sig 2: %x %x %x", magic[0],magic[1],magic[2]);

  // NETCDF
  if (magic[0]==0x43 && magic[1]==0x44 && magic[2]==0x46) {
    if (debug>0) mprintf("  NETCDF file\n");
    if (compressType!=NO_COMPRESSION) {
      mprintf("Error: Compressed NETCDF files are not currently supported.\n");
      return 1;
    }
    // Determine whether this is a trajectory or restart Netcdf from the Conventions
    CheckConventions = GetNetcdfConventions(filename);
    if (CheckConventions==NULL) return 1;
    if (strcmp(CheckConventions, "AMBER")==0) 
      fileFormat=AMBERNETCDF;
    else if (strcmp(CheckConventions,"AMBERRESTART")==0)
      fileFormat=AMBERRESTARTNC;
    else {
      mprintf("Error: Netcdf File %s: Unrecognized conventions \"%s\".\n",
              filename,CheckConventions);
      mprintf("       Expected \"AMBER\" or \"AMBERRESTART\".\n");
      fileFormat=UNKNOWN_FORMAT;
      free(CheckConventions);
      return 1;
    }
    free(CheckConventions);
    return 0;
  }

  // ---------- ID by file characteristics; read the first two lines ----------
  // Initialize buffers to NULL
  memset(buffer1,' ',BUFFER_SIZE);
  memset(buffer2,' ',BUFFER_SIZE);
  buffer1[0]='\0';
  buffer2[0]='\0';
  IO->Open(filename,"rb"); // NOTE: Err Check
  IO->Gets(buffer1,BUFFER_SIZE);
  IO->Gets(buffer2,BUFFER_SIZE);
  IO->Close();

  // Check for terminal CR before newline, indicates DOS file
  i = strlen(buffer1);
  if ( i>1 ) {
    if (buffer1[ i - 2 ] == '\r') {
      if (debug>0) mprintf("  [DOS]");
      isDos=1;
    }
  }

  // If both lines have PDB keywords, assume PDB
  if (isPDBkeyword(buffer1) && isPDBkeyword(buffer2)) {
    if (debug>0) mprintf("  PDB file\n");
    fileFormat=PDBFILE;
    return 0;
  }

  // If either buffer contains a TRIPOS keyword assume Mol2
  // NOTE: This will fail on tripos files with extensive header comments.
  //       A more expensive check for mol2 files is below.
  if (strncmp(buffer1,"@<TRIPOS>", 9)==0 ||
      strncmp(buffer2,"@<TRIPOS>", 9)==0) {
    if (debug>0) mprintf("  TRIPOS MOL2 file\n");
    fileFormat=MOL2FILE;
    return 0;
  }

  // If the %VERSION and %FLAG identifiers are present, assume amber parm
  if (strncmp(buffer1,"%VERSION",8)==0 && strncmp(buffer2,"%FLAG",5)==0) {
    if (debug>0) mprintf("  AMBER TOPOLOGY file\n");
    fileFormat=AMBERPARM;
    return 0;
  }

  // If the first 3 chars are P S F, assume char PSF
  if (strncmp(buffer1,"PSF",3)==0) {
    if (debug>0) mprintf("  CHARMM PSF file\n");
    fileFormat=CHARMMPSF;
    return 0;
  }

  // If the second 5 chars are C O R D, assume charmm DCD
  if (strncmp(buffer1+4,"CORD",4)==0) {
    if (debug>0) mprintf("  CHARMM DCD file\n");
    fileFormat=CHARMMDCD;
    return 0;
  }

  // Amber Restart
  // Check for an integer (I5) followed by 0-2 scientific floats (E15.7)
  if (strlen(buffer2)<=36) {
    //mprintf("DEBUG: Checking restart.\n");
    //mprintf("DEBUG: buffer2=[%s]\n",buffer2);
    for (i=0; i<5; i++) {
      if (!isspace(buffer2[i]) && !isdigit(buffer2[i])) break;
      //mprintf("DEBUG:    %c is a digit/space.\n",buffer2[i]);
    }
    //mprintf("DEBUG: i=%i\n");
    //if ( i==5 && strchr(buffer2,'E')!=NULL ) {
    if ( i==5 ) {
      if (debug>0) mprintf("  AMBER RESTART file\n");
      fileFormat=AMBERRESTART;
      return 0;
    }
  }

  // If first line is 81 bytes and the second line has 12 numbers in
  // 12I6 format, assume old-style Amber topology
  // NOTE: Could also be less than 81? Only look for 12 numbers?
  if ((int)strlen(buffer1)==81+isDos) {
    if ( sscanf(buffer2,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i", iamber, iamber+1,
                iamber+2, iamber+3, iamber+4, iamber+5, iamber+6, 
                iamber+7, iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
    {
      if (debug>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
      fileFormat=OLDAMBERPARM;
      return 0;
    }
  }

  // Check if second line contains REMD/HREMD, Amber Traj with REMD header
  if ( strncmp(buffer2,"REMD",4)==0 ||
       strncmp(buffer2,"HREMD",5)==0   ) 
  {
    if (debug>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
    fileFormat=AMBERTRAJ;
    return 0;
  }

  // Check if we can read at least 3 coords of width 8, Amber trajectory
  if ( sscanf(buffer2, "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 ) {
    if (debug>0) mprintf("  AMBER TRAJECTORY file\n");
    fileFormat=AMBERTRAJ;
    return 0;
  }

  // OLD STYLE CHECK: Would not recognize 1 atom trajectories, dependent on line size.
  // If the second line is 81 bytes, assume Amber Traj
  // If second line is 42 bytes, check for Amber Traj with REMD
  // NOTE: Check for digit every 8 char?
  // NOTE: This wont work for trajectories with <3 coords.
/*  if ((int)strlen(buffer2)==81+isDos) {
    //if ( sscanf(buffer2, "%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2,
    //            TrajCoord+3, TrajCoord+4, TrajCoord+5, TrajCoord+6, TrajCoord+7, TrajCoord+8,
    //            TrajCoord+9) == 10 )
    // Make sure we can read at least 3 coords of width 8
    
  } else if (strlen(buffer2)==42) {
    
  }
*/

  // ---------- MORE EXPENSIVE CHECKS ----------
  // Reopen and scan for Tripos mol2 molecule section
  // 0 indicates section found.
  IO->Open(filename,"r");
  if (!Mol2ScanTo(this, MOLECULE)) {
    if (debug>0) mprintf("  TRIPOS MOL2 file\n");
    fileFormat=MOL2FILE;
    IO->Close();
    return 0;
  }
  IO->Close();

  // ---------- EXPERIMENTAL ----------
  // If the file format is still undetermined and the file name is conflib.dat,
  // assume this is a conflib.dat file from LMOD. Cant think of a better way to
  // detect this since there is no magic number but the file is binary.
  if ( strcmp(basefilename,"conflib.dat")==0 ) {
    mprintf("  LMOD CONFLIB file\n");
    fileFormat=CONFLIB;
    return 0;
  }

  // Unidentified file
  mprintf("  Warning: %s: UNKNOWN FILE FORMAT.\n",filename);
  return 1; 
}

