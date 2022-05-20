// FileIO_Gzip: Gzip file operations
#ifdef HASGZ
#include <cstdio>
#include "FileIO_Gzip.h" // FileIO.h, zlib.h
#include "CpptrajStdio.h"

// CONSTRUCTOR
FileIO_Gzip::FileIO_Gzip() : fp_(NULL) {}

// DESTRUCTOR
FileIO_Gzip::~FileIO_Gzip() {
  Close();
}

//const unsigned int FileIO_Gzip::GZ_BUF_SIZE = 65536;
//const unsigned int FileIO_Gzip::GZ_BUF_SIZE = 131072;

// FileIO_Gzip::Open()
int FileIO_Gzip::Open(const char *filename, const char *mode) {
  if (filename == 0) return 1;
  fp_ = gzopen(filename, mode);
  if (fp_ == NULL) return 1;
  // Set gzip buffer size
  //if ( gzbuffer(fp_, GZ_BUF_SIZE)!=0 ) return 1;
  return 0;
}

// FileIO_Gzip::Close()
int FileIO_Gzip::Close() {
  if (fp_ != NULL)
    gzclose(fp_);
  fp_ = NULL;
  return 0;
}

// FileIO_Gzip::Size()
/** Gzip files include the uncompressed size % 2^32 in the last 4 bytes 
  * of the file. Return this value. This is used for example when attempting
  * to determine the number of frames in a gzip compressed amber traj.
  */
off_t FileIO_Gzip::Size(const char *filename) {
  FILE *infile;
  unsigned char b1,b2,b3,b4;
  off_t val,temp;

  if (filename==NULL) return -1;
  if ( (infile = fopen(filename,"rb"))==NULL ) {
    mprintf("Error: FileIO_Gzip::Size: Could not open %s for reading.\n",filename);
    return -1L;
  }

  // Place 4 bytes from the end
  fseek(infile, -4, SEEK_END);

  b1=0; b2=0; b3=0; b4=0;
  if (fread(&b4,1,1,infile) != 1) return -1L;
  if (fread(&b3,1,1,infile) != 1) return -1L;
  if (fread(&b2,1,1,infile) != 1) return -1L;
  if (fread(&b1,1,1,infile) != 1) return -1L;

  val = 0;
  temp = (off_t) b1;
  temp <<= 24;
  val = val | temp;
  temp = (off_t) b2;
  temp <<= 16;
  val = val | temp;
  temp = (off_t) b3;
  temp <<= 8;
  val = val | temp;
  temp = (off_t) b4;
  val = val | temp;

  fclose(infile);

  //fprintf(stdout,"FileIO_Gzip::Size: Uncompressed size of %s: %lu\n",filename,val);

  return val;
}

// FileIO_Gzip::Read()
// NOTE: gzread returns 0 on EOF, -1 on error
int FileIO_Gzip::Read(void *buffer, size_t num_bytes) {
  int n_bytes_read = gzread(fp_, buffer, num_bytes);
  if (n_bytes_read < 1) {
    int gzerr = 0;
    const char* gzerrmsg = gzerror(fp_, &gzerr);
    mprinterr("Internal Error: FileIO_Gzip::Read: %s\n", gzerrmsg);
  }
  return n_bytes_read;
}

// FileIO_Gzip::Write()
int FileIO_Gzip::Write(const void *buffer, size_t num_bytes) {
  int n_bytes_written = gzwrite(fp_, buffer, num_bytes);
  if ((size_t)n_bytes_written != num_bytes) {
    int gzerr = 0;
    const char* gzerrmsg = gzerror(fp_, &gzerr);
    mprinterr("Internal Error: FileIO_Gzip::Write: %s\n", gzerrmsg);
    return 1;
  }
  return 0;
}

// FileIO_Gzip::Seek()
int FileIO_Gzip::Seek(off_t offset) {
  z_off_t zipOffset;
 
  //if (origin == SEEK_END) return 1; 
  zipOffset = (z_off_t)offset;
  z_off_t uncompressed_pos = gzseek(fp_, zipOffset, SEEK_SET) < 0;
  if (uncompressed_pos < 0) {
    mprinterr("Internal Error: FileIO_Gzip::Seek offset = %li.\n", (long int)offset);
    return 1;
  }
  return 0;
}

// FileIO_Gzip::Rewind()
int FileIO_Gzip::Rewind() {
  gzrewind(fp_);
  return 0;
}

// FileIO_Gzip::Tell()
off_t FileIO_Gzip::Tell() {
  z_off_t zipOffset;
  
  zipOffset = gztell(fp_);
  return (off_t) zipOffset;
}

// FileIO_Gzip::Gets()
int FileIO_Gzip::Gets(char *str, int num) {
  if ( gzgets(fp_,str,num) == NULL )
    return 1;
  else
    return 0;
}
#endif
