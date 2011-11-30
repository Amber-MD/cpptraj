// FileIO_Bzip2: Bzip2 file operations
#ifdef HASBZ2
#include <cstring>
#include <cstdlib>
#include "FileIO_Bzip2.h" // FileIO.h, cstdio, bzlib.h
#include "CpptrajStdio.h"

// CONSTRUCTOR
FileIO_Bzip2::FileIO_Bzip2() {
  //fprintf(stderr,"FileIO_Bzip2 CONSTRUCTOR\n");
  fp = NULL;
  infile = NULL;
  err = BZ_OK;
  isBzread=true;
  bzfilename=NULL;
  bzmode=NULL;
  position=0L;
}

// DESTRUCTOR
FileIO_Bzip2::~FileIO_Bzip2() {
  //fprintf(stderr,"FileIO_Bzip2 DESTRUCTOR\n");
  if (fp!=NULL || infile!=NULL) this->Close();
  if (bzfilename!=NULL) free(bzfilename);
  if (bzmode!=NULL) free(bzmode);
}

// FileIO_Bzip2::BZerror()
/** Return a string corresponding to the current value of err.
  */
const char *FileIO_Bzip2::BZerror() {
  switch (err) {
    case BZ_OK : return "BZ_OK";
    case BZ_PARAM_ERROR : return "BZ_PARAM_ERROR";
    case BZ_SEQUENCE_ERROR : return "BZ_SEQUENCE_ERROR";
    case BZ_IO_ERROR: return "BZ_IO_ERROR";
    case BZ_UNEXPECTED_EOF : return "BZ_UNEXPECTED_EOF";
    case BZ_DATA_ERROR : return "BZ_DATA_ERROR";
    case BZ_DATA_ERROR_MAGIC : return "BZ_DATA_ERROR_MAGIC";
    case BZ_MEM_ERROR : return "BZ_MEM_ERROR";
    case BZ_STREAM_END : return "BZ_MEM_ERROR";
  }
  return "Unknown Bzip2 error";
}

// FileIO_Bzip2::Open()
/** Open the given file as a bzip2 file. The mode and filename are stored
  * in case rewind is called (Bzip2 routines do not have a rewind so the
  * file must be closed and reopened).
  */
// NOTES from Bzip2 docs:
// BZFILE *BZ2_bzReadOpen( int *bzerror, FILE *f, int verbosity, int small,
//                         void *unused, int nUnused );
//   If small is 1, the library will try to decompress using less memory, at 
//   the expense of speed. BZ2_bzRead will decompress the nUnused bytes starting 
//   at unused, before starting to read from the file f. At most BZ_MAX_UNUSED 
//   bytes may be supplied like this. If this facility is not required, you 
//   should pass NULL and 0 for unused and nUnused respectively.
// BZFILE *BZ2_bzWriteOpen( int *bzerror, FILE *f, int blockSize100k, 
//                          int verbosity, int workFactor );
//   Parameter blockSize100k specifies the block size to be used for compression.
//   It should be a value between 1 and 9 inclusive, and the actual block size 
//   used is 100000 x this figure. 9 gives the best compression but takes most 
//   memory. Parameter verbosity should be set to a number between 0 and 4 
//   inclusive. 0 is silent. Parameter workFactor controls how the compression 
//   phase behaves when presented with worst case, highly repetitive, input data.
//   If compression runs into difficulties caused by repetitive data, the library 
//   switches from the standard sorting algorithm to a fallback algorithm. The 
//   fallback is slower than the standard algorithm by perhaps a factor of three, 
//   but always behaves reasonably, no matter how bad the input. Lower values of 
//   workFactor reduce the amount of effort the standard algorithm will expend 
//   before resorting to the fallback. You should set this parameter carefully; 
//   too low, and many inputs will be handled by the fallback algorithm and so 
//   compress rather slowly, too high, and your average-to-worst case compression
//   times can become very large. The default value of 30 gives reasonable 
//   behaviour over a wide range of circumstances.
int FileIO_Bzip2::Open(const char *filename, const char *mode) {
  // Store filename and mode - reallocate in case of reopen
  if (bzfilename!=filename) {
    bzfilename = (char*) realloc(bzfilename, (strlen(filename)+1) * sizeof(char));
    strcpy(bzfilename, filename);
  }
  if (bzmode!=mode) {
    bzmode     = (char*) realloc(bzmode,     (strlen(mode)+1    ) * sizeof(char));
    strcpy(bzmode, mode);
  }

  // DEBUG
  //mprintf("DEBUG: FileIO_Bzip2::Open(%s,%s)\n",filename,mode);

  fp = fopen(filename, mode);
  if (fp==NULL) {
    mprintf("Error: FileIO_Bzip2::Open: Could not open %s with mode %s\n",filename,mode);
    return 1;
  }

  switch ( mode[0] ) {
    case 'r' : 
      //mprintf("DEBUG: Calling bzReadOpen\n");
      infile = BZ2_bzReadOpen( &err, fp, 1, 0, NULL, 0); 
      isBzread=true;  
      break;
    case 'w' : 
      //mprintf("DEBUG: Calling bzWriteOpen\n");
      infile = BZ2_bzWriteOpen( &err, fp, 9, 0, 30);
      isBzread=false; 
      break;
    case 'a' : 
      mprintf("Error: FileIO_Bzip2::Open: Append not supported for Bzip2.\n");
      return 1; // No append for Bzip2
    default: return 1; 
  }

  if (err != BZ_OK) {
    mprintf("Error: FileIO_Bzip2::Open: Could not BZOPEN %s with mode %s\n",filename,mode);
    return 1;
  }

  if (infile==NULL) return 1;
  //mprintf("DEBUG: BZIP2 Opened %s with mode %s\n",filename,mode);
  position=0L;
  return 0;
}

// FileIO_Bzip2::Close()
int FileIO_Bzip2::Close() {
  if (infile!=NULL) {
    if (isBzread) {
      //mprintf("DEBUG: BZ2_bzReadClose\n");
      BZ2_bzReadClose(&err, infile);
    } else {
      //mprintf("DEBUG: BZ2_bzWriteClose\n");
      BZ2_bzWriteClose(&err, infile, 0, NULL, NULL);
    }
    infile=NULL;
  }
  
  if (fp!=NULL) fclose(fp);
  fp=NULL;
  return 0;
}

// FileIO_Bzip2::Size()
/** Since the uncompressed size of Bzip files is not stored anywhere in the file
  * need to read every possible byte in the file, which can be VERY slow. 
  * NOTE: The input filename is currently IGNORED.
  * NOTE: This can be ridiculously time consuming for large bzip files, so
  *       just return 0. 
  */
#define BUFINSIZE 10240
off_t FileIO_Bzip2::Size(char *filename) {
  //off_t fileSize, numread;
  //char bufIn[BUFINSIZE];
//  char Scan;

  if (filename==NULL) return -1L;
  //fileSize=0L;
  return 0L;
/*
  // Check that the file being checked is the currently open file.
  // NOTE: Unnecessary?
  //if (strcmp(filename, bzfilename)!=0) {
  //  mprintf("ERROR: FileIO_Bzip2::Size: Checking file %s, open file is %s!\n",
  //          filename,bzfilename);
  //  return -1L;
  //}

  // Open the file
  if (infile==NULL) {
    mprintf("FileIO_Bzip2::Size: Opening %s\n",filename);
    if (this->Open(filename,"rb")) return -1L;
  }

  // Read all chars in file 10K bytes at a time.
  // NOTE: Use sizeof(char)??
  // NOTE: Larger buffer? Dynamically allocate?
//  bufIn = (char*) malloc(10240 * sizeof(char));
  while ( (numread = (off_t) this->Read(bufIn, 1, BUFINSIZE))!=-1 )
    fileSize += numread;
//  free(bufIn);
//  while ( this->Read(&Scan, 1, 1)!=-1 )
//    fileSize = fileSize + 1L;

  // Close file
  this->Close();

  //mprintf("FileIO_Bzip2::Size: Uncompressed size of %s: %lu\n",filename,fileSize);

  return fileSize;
*/
}
#undef BUFINSIZE

// FileIO_Bzip2::Read()
/** Read size*count bytes from bzip2file stream. Return number of bytes read.
  * If an error occurs or no more bytes to be read return -1;
  * Dont attempt to read if error bit is set.
  */
int FileIO_Bzip2::Read(void *buffer, size_t size, size_t count) {
  //size_t numread;
  int numread;
  int expectedread;
  
  //if (err!=BZ_OK) return -1;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  mprintf("Error: FileIO_Bzip2::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  expectedread = (int) size;
  expectedread *= (int) count;
  numread = BZ2_bzRead(&err, infile, buffer, expectedread);

  // Update position
  position = position + ((off_t) numread);

  if (numread != expectedread) {
    if (err!=BZ_OK && err!=BZ_STREAM_END) {
      mprintf( "Error: FileIO_Bzip2::Read: BZ2_bzRead error: [%s]\n",this->BZerror());
      mprintf( "                        size=%lu  count=%lu\n",size,count);
    }
    return -1;
  }
  //mprintf( "DEBUG: After FileIO_Bzip2::Read: [%s] position %li\n",this->BZerror(),position);

  return numread;
}

// FileIO_Bzip2::Write()
int FileIO_Bzip2::Write(void *buffer, size_t size, size_t count) {
  //size_t numwrite;
  int numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  mprintf("Error: FileIO_Bzip2::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  numwrite = size * count;
  BZ2_bzWrite ( &err, infile, buffer, numwrite );

  // Update position
  position = position + ((off_t)numwrite);

  if (err == BZ_IO_ERROR) { 
    mprintf( "Error: FileIO_Bzip2::Write: BZ2_bzWrite error\n");
    return 1;
  }

  return 0;
}

// FileIO_Bzip2::Seek
/** Since a true seek is not really possible with bzip2, scan 1 char at 
  * a time until desired position achieved.
  */
// NOTE: Scan in blocks?
int FileIO_Bzip2::Seek(off_t offset) {
  off_t seekTo;
  char Scan;
  // Determine place to seek to
  //switch (origin) {
  //  case SEEK_SET : seekTo = offset; break;
  //  case SEEK_CUR : seekTo = position + offset; break;
  //  case SEEK_END : 
  //    mprintf("Error: FileIO_Bzip2::Seek: Seek to END not supported (%s).\n",bzfilename);
  //  default : return 1;
  //}
  seekTo = offset;

  //mprintf("DEBUG: FileIO_Bzip2::Seek: %s %li -> %li, ",bzfilename,position,seekTo);

  // If place to seek to is earlier than current position need to reopen
  if (seekTo<position) 
    this->Rewind();

  // Read chars until position achieved
  while (position<seekTo) {
    if (this->Read(&Scan,1,1)==-1) break;
  }

  //mprintf("%li\n",position);

  return 0;
}

// FileIO_Bzip2::Rewind()
/** Close and reopen.
  */
int FileIO_Bzip2::Rewind() {
  if (bzfilename==NULL || bzmode==NULL) return 1;
  this->Close();
  this->Open(bzfilename,bzmode);
  return 0;
}

// FileIO_Bzip2::Tell()
// NOTE: Tell not possible with bzip2. Use position.
off_t FileIO_Bzip2::Tell() {
  return position;
}

// FileIO_Bzip2::Gets()
/** Analogous to fgets, reads characters from stream and stores them as a C 
  * string into str until (num-1) characters have been read or either a newline
  * or the End-of-File is reached, whichever comes first.
  * A newline character makes fgets stop reading, but it is considered a valid 
  * character and therefore it is included in the string copied to str.
  * A null character is automatically appended in str after the characters read
  * to signal the end of the C string.
  */
int FileIO_Bzip2::Gets(char *str, int num) {
  int i;
  //mprintf("DEBUG: FileIO_Bzip2::Gets: num=%i\n",num);
  // Try to read num chars. If newline encountered exit
  if (num<=1) return 1;
  i=0;
  while ( this->Read(str+i, 1, 1)!=-1 ) {
    i++;
    if (i==num-1) break;
    if (str[i-1]=='\n') break;
  }
  // If nothing read return 1
  if (i==0) return 1;
  // i should be at num or 1 after newline; append NULL char
  str[i] = '\0';
  //mprintf("DEBUG: FileIO_Bzip2::Gets: num=%i i=%i [%s]\n",num,i,str);
  //mprintf( "DEBUG: After FileIO_Bzip2::Gets: position %li\n",position);
  return 0;
}
#endif
