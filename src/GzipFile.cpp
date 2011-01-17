// GzipFile: Gzip file operations
#ifdef HASGZ
#include <cstdio>
#include "GzipFile.h" // BaseFileIO.h, zlib.h
#include "CpptrajStdio.h"

// CONSTRUCTOR
GzipFile::GzipFile() {
  fp = NULL;
}

// DESTRUCTOR
GzipFile::~GzipFile() {
  if (fp!=NULL) this->Close();
}

/*
 * GzipFile::Open()
 */
int GzipFile::Open(const char *filename, const char *mode) {
  fp = gzopen(filename, mode);
  if (fp==NULL) return 1;
  return 0;
}

/*
 * GzipFile::Close()
 */
int GzipFile::Close() {
  if (fp!=NULL) gzclose(fp);
  fp=NULL;
  return 0;
}

/*
 * GzipFile::Size()
 */
long long int GzipFile::Size(char *filename) {
  FILE *infile;
  unsigned char b1,b2,b3,b4;
  long long int val,temp;

  if (filename==NULL) return -1;
  if ( (infile = fopen(filename,"rb"))==NULL ) {
    mprintf("Error: GzipFile::Size: Could not open %s for reading.\n",filename);
    return -1L;
  }

  // Place 4 bytes from the end
  fseek(infile, -4, SEEK_END);

  b1=0; b2=0; b3=0; b4=0;
  fread(&b4,1,1,infile);
  fread(&b3,1,1,infile);
  fread(&b2,1,1,infile);
  fread(&b1,1,1,infile);

  val = 0;
  temp = (long long int) b1;
  temp <<= 24;
  val = val | temp;
  temp = (long long int) b2;
  temp <<= 16;
  val = val | temp;
  temp = (long long int) b3;
  temp <<= 8;
  val = val | temp;
  temp = (long long int) b4;
  val = val | temp;

  fclose(infile);

  //fprintf(stdout,"GzipFile::Size: Uncompressed size of %s: %lli\n",filename,val);

  return val;
}

/*
 * GzipFile::Read()
 */
int GzipFile::Read(void *buffer, size_t size, size_t count) {
  //size_t numread;
  int numread;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: GzipFile::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  numread = gzread(fp, buffer, size * count);
  if (numread == -1) return -1;

  // NOTE: Check for errors here.
  return numread;
}

/*
 * GzipFile::Write()
 */
int GzipFile::Write(void *buffer, size_t size, size_t count) {
  //size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: GzipFile::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  if ( gzwrite(fp, buffer, size * count)==0 ) return 1;
  // NOTE: Check for errors here.
  return 0;
}

/*
 * GzipFile::Seek()
 */
int GzipFile::Seek(off_t offset) {
  z_off_t zipOffset;
 
  //if (origin == SEEK_END) return 1; 
  zipOffset=(z_off_t) offset;
  if ( gzseek(fp, zipOffset, SEEK_SET) < 0) return 1;
  return 0;
}

/*
 * GzipFile::Rewind()
 */
int GzipFile::Rewind() {
  gzrewind(fp);
  return 0;
}

/*
 * GzipFile::Tell()
 */
off_t GzipFile::Tell() {
  z_off_t zipOffset;
  
  zipOffset = gztell(fp);
  return (off_t) zipOffset;
}

/*
 * GzipFile::Gets()
 */
int GzipFile::Gets(char *str, int num) {
  if ( gzgets(fp,str,num) == NULL )
    return 1;
  else
    return 0;
}
#endif
