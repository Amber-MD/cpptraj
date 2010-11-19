// PtrajFile 
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <glob.h> // For tilde expansion
#include "PtrajFile.h"
// File Types
#include "StdFile.h"
#ifdef HASGZ
#  include "GzipFile.h"
#endif
#include "MpiFile.h"
#ifdef HASBZ2
#  include "Bzip2File.h"
#endif

//typedef char enumToken[30];
const PtrajFile::enumToken PtrajFile::FileFormatList[9] = {
  "UNKNOWN_FORMAT", "PDBFILE", "AMBERTRAJ", "AMBERNETCDF", "AMBERPARM", 
  "DATAFILE", "AMBERRESTART", "AMBERREMD", "XMGRACE"
};
const PtrajFile::enumToken PtrajFile::FileTypeList[6] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};

// CONSTRUCTOR
PtrajFile::PtrajFile() {
  IO=NULL;
  isOpen=0;
  uncompressed_size=0UL;
  compressType=NONE;
  access=READ;
  debug=0;
  fileFormat=UNKNOWN_FORMAT;
  fileType=STANDARD;
  filename=NULL;
  basefilename=NULL;
  Ext=NULL;
  isDos=0;
}

// Return enumerated type for FileFormat
char *PtrajFile::Format() {
  return (char*)FileFormatList[fileFormat];
}

// Return enumerated type for FileType
char *PtrajFile::Type() {
  return (char*)FileTypeList[fileType];
}

// DESTRUCTOR
PtrajFile::~PtrajFile() {
   //fprintf(stderr,"PTRAJFILE DESTRUCTOR\n");
   CloseFile();
   delete IO;
   if (filename!=NULL) free(filename);
   if (basefilename!=NULL) free(basefilename);
   if (Ext!=NULL) free(Ext);
}

/*
 * isPDBkeyword()
 * Return 1 if the first 6 chars of buffer match a PDB keyword
 * NOTE: Should really be in PDBFileRoutines, but put here since
 * its currently only used by the PtrajFile class.
 */
int isPDBkeyword(char *buffer) {
  if (strncmp(buffer,"HEADER",6)==0) return 1;
  if (strncmp(buffer,"TITLE ",6)==0) return 1;
  if (strncmp(buffer,"COMPND",6)==0) return 1;
  if (strncmp(buffer,"ATOM  ",6)==0) return 1;
  if (strncmp(buffer,"CRYST1",6)==0) return 1;
  if (strncmp(buffer,"REMARK",6)==0) return 1;
  if (strncmp(buffer,"MODEL ",6)==0) return 1;
  return 0;
}

/* PtrajFile::OpenFile()
 * Open the file. If already open, reopen.
 */
int PtrajFile::OpenFile() {
  if (isOpen) CloseFile();

  switch (access) {
    case READ:
      if ( IO->Open(filename, "rb")  ) { // NOTE: use rb as mode instead?
        fprintf(stdout,"Could not open %s for reading.\n",filename);
        return 1;
      }
      if (debug>0) fprintf(stdout,"Opened %s for reading.\n",filename);
      break;
    case APPEND:
      if ( IO->Open(filename, "ab") ) {
        fprintf(stdout,"Could not open %s for appending.\n",filename);
        return 1;
      }
      if (debug>0) fprintf(stdout,"Opened %s for appending.\n",filename);
      break;
    case WRITE:
      if ( IO->Open(filename, "wb") ) { // NOTE: Use wb as mode?
        fprintf(stdout,"Could not open %s for writing.\n",filename);
        return 1;
      }
      if (debug>0) fprintf(stdout,"Opened %s for writing.\n",filename);
      break;
  }
      
  isOpen=1;
  return 0;
}

/* PtrajFile::CloseFile()
 * Close the file.
 */
void PtrajFile::CloseFile() {
  if (isOpen) {
    IO->Close();
    isOpen=0;
  }
}

/* PtrajFile::determineType()
 * If the file type is unknown attempt to determine it from filename extension.
 * Default to standard.
 */
void PtrajFile::determineType() {

  if (fileType!=UNKNOWN_TYPE) return;

  // No extension
  if (Ext==NULL) {
    fileType=STANDARD;
    return;
  }

  if      ( strcmp(Ext,".gz")==0  ) fileType=GZIPFILE;
  else if ( strcmp(Ext,".bz2")==0 ) fileType=BZIP2FILE;
  else fileType=STANDARD;
}

/*
 * PtrajFile::determineFormat()
 * If the file format is unknown attempt to determine from filename extension.
 * Default to datafile.
 */
void PtrajFile::determineFormat() {

  if (fileFormat!=UNKNOWN_FORMAT) return;

  // No extension
  if (Ext==NULL) return;

  if      ( strcmp(Ext,".dat")==0 ) fileFormat=DATAFILE;
  else if ( strcmp(Ext,".agr")==0 ) fileFormat=XMGRACE;
  else fileFormat=DATAFILE;
}  

/*
 * PtrajFile::SetBaseFilename()
 * Strip leading path from input filename. Use strtok routine to separate 
 * filename by / and use the last string as the base filename. Internal 
 * filename is not used since strtok modifies the char array.
 * Also determine the file extension.
 */
void PtrajFile::SetBaseFilename() {
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
    printf("DEBUG: SetBaseFilename: Filename is %s\n",filename);
    printf("                        Base filename is %s\n",basefilename);
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
      fprintf(stdout,"PtrajFile: No extension.\n");
    return;
  } 
  Ext = (char*) malloc( (strlen(basefilename+i) + 1) * sizeof(char));
  strcpy(Ext, basefilename+i);
  if (debug>0)
    fprintf(stdout,"PtrajFile: Extension= %s  Length= %lu\n",Ext,strlen(basefilename+i));
}

/* 
 * PtrajFile::SetupFile()
 * Sets the file name, access type (R/W/A), file type and file format (for 
 * WRITE), and debug level. If called with READ or append file type and format
 * will be determined by SetupRead. If called with WRITE the given type and 
 * format will be used; if called with UNKNOWN type and format they will
 * try to be determined by the given file extension. 
 * Can be called with NULL for write, this will write to stdout.
 */
int PtrajFile::SetupFile(char *filenameIn, AccessType accessIn, 
                         FileFormat fileFormatIn, FileType fileTypeIn, 
                         int debugIn) {
  glob_t globbuf;

  // DEBUG
  if (debug>1) {
    printf("PTRAJFILE::SETUPFILE0: %s called with format %s and type %s\n",filenameIn,
           FileFormatList[fileFormatIn],FileTypeList[fileTypeIn]);
  }

  debug=debugIn;
  // Store filename
  if (filenameIn!=NULL) {
    filename=(char*) malloc( (strlen(filenameIn)+1) * sizeof(char) );
    strcpy(filename,filenameIn);
  }
  // On read or append do tilde expansion and store new filename
  if (accessIn!=WRITE) {
    // If no filename this is an error.
    if (filename==NULL) {
      fprintf(stdout,"Error: PtrajFile::SetupFile: NULL filename specified on READ or APPEND\n");
      return 1;
    }
    globbuf.gl_offs = 1;
    if ( glob(filename, GLOB_TILDE, NULL, &globbuf)!=0 ) return 1; 
    if (debug>1) fprintf(stdout,"  GLOB(0): [%s]\n",globbuf.gl_pathv[0]);
    filename=(char*) realloc( filename, (strlen(globbuf.gl_pathv[0])+1) * sizeof(char));
    strcpy(filename, globbuf.gl_pathv[0]);
    globfree(&globbuf);
  }
  // Store base filename and determine filename extension
  this->SetBaseFilename();

  access=accessIn;
  switch (access) {
    case APPEND:
    case READ:  
      // SetupRead determines type and format
      if ( SetupRead()  ) return 1; 
      break;
    case WRITE: 
      // If file type is not specified, try to determine from filename
      // NOTE: Also determine format from name??
      fileType=fileTypeIn;
      fileFormat=fileFormatIn;
      this->determineType();
      this->determineFormat();
      if ( SetupWrite() ) return 1; 
      break;
    default: return 1;
  }
  if (debug>0)
    printf("PtrajFile::SetupFile: %s is format %s and type %s\n",filename,
           FileFormatList[fileFormat],FileTypeList[fileType]);
  return 0;
}

/* 
 * PtrajFile::SetupWrite()
 * Set up file with specified type for writing. fileFormat is set by SetupFile.
 * NOTE: Determine if compression is requested, either by arg or from name
 */
int PtrajFile::SetupWrite() {

  //fprintf(stdout,"DEBUG: Setting up write file %s with format %s\n",filename,fileformatIn);
  // Eventually allow other file types
  //fileType=STANDARD;
  switch (fileType) {
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new GzipFile(); 
#else
      fprintf(stdout,"Error: SetupWrite(%s):\n",filename);
      fprintf(stdout,"       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      IO = new Bzip2File();
#else
      fprintf(stdout,"Error: SetupWrite(%s):\n",filename);
      fprintf(stdout,"       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#endif
    break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    case STANDARD  : IO = new StdFile();  break;
    case MPIFILE   : IO = new MpiFile();  break;
    default : 
      fprintf(stdout,"PtrajFile::SetupWrite: Unrecognized file type.\n");
      return 1;
      break;
  }

  return 0;
}

/*
 * PtrajFile::SetupRead() 
 * Open the file specified by filename for READ or APPEND access. Attempt to 
 * identify the file type and format.
 */
int PtrajFile::SetupRead() {
  unsigned char magic[3];
  char buffer1[BUFFER_SIZE], buffer2[BUFFER_SIZE];
  float TrajCoord[10];
  int i;

  //fprintf(stdout,"DEBUG: Setting up read file %s\n",filename);
  // Get basic file information
  // An error here means file probably doesnt exist. Dont print an error at 
  // basic debug level since this could also be used to test if file exists.
  if (stat(filename, &frame_stat) == -1) {
    if (debug>0) {
      fprintf(stdout, "ERROR: PtrajFile::SetupRead: Could not find file status for %s\n", filename);
      perror("     Error from stat: ");
    }
    return 1;
  }

  // Start off every file as a standard file - may need to change for MPI
  IO = new StdFile();

  // ID by magic number - open for binary read access
  if ( IO->Open(filename, "rb") ) { 
    fprintf(stdout,"Could not open %s for hex signature read.\n",filename);
    return 1;
  }

  // Read first 3 bytes
  memset(magic,0,3*sizeof(unsigned char));
  IO->Read(magic  ,1,1);
  IO->Read(magic+1,1,1);
  IO->Read(magic+2,1,1);
  IO->Close();
  if (debug>0) fprintf(stdout,"File: %s: Hex sig: %x %x %x", filename,
                       magic[0],magic[1],magic[2]);

  // Check compression
  if ((magic[0]==0x1f) && (magic[1]==0x8b) && (magic[2]==0x8)) {
    if (debug>0) fprintf(stdout,", Gzip file.\n");
    compressType=GZIP;
    fileType=GZIPFILE;
  } else if ((magic[0]==0x42) && (magic[1]==0x5a) && (magic[2]==0x68)) {
    if (debug>0) fprintf(stdout,", Bzip2 file.\n");
    compressType=BZIP2;
    fileType=BZIP2FILE;
  } else if ((magic[0]==0x50) && (magic[1]==0x4b) && (magic[2]==0x3)) {
    if (debug>0) fprintf(stdout,", Zip file.\n");
    compressType=ZIP;
    fileType=ZIPFILE;
  } else {
    if (debug>0) fprintf(stdout,", No compression.\n");
  }

  // Appending and compression not supported.
  if (access==APPEND && compressType!=NONE) {
    fprintf(stdout,"Error: Appending to compressed files is not supported.\n");
    return 1;
  }

  // Assign the appropriate IO type based on the file
  delete (StdFile*) IO;
  IO = NULL;
  switch (fileType) {
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new GzipFile(); 
#else
      fprintf(stdout,"Error: SetupRead(%s):\n",filename);
      fprintf(stdout,"       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE : 
#ifdef HASBZ2
      IO = new Bzip2File(); 
#else
      fprintf(stdout,"Error: SetupRead(%s):\n",filename);
      fprintf(stdout,"       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#endif
      break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    default        : IO = new StdFile();  break;
  }

  // If the file is compressed, get the uncompressed size
  // For standard files this just returns 0UL
  // Standard file size is in the frame_stat struct
  uncompressed_size = IO->Size(filename);

  // Determine format
  // Read first 3 bytes again to determine format by magic number
  IO->Open(filename,"rb"); // NOTE: Err Check
  memset(magic,0,3*sizeof(unsigned char));
  IO->Read(magic  ,1,1);
  IO->Read(magic+1,1,1);
  IO->Read(magic+2,1,1);
  IO->Close();
  if (debug>0) fprintf(stdout,"File: %s: Hex sig 2: %x %x %x", filename,
                       magic[0],magic[1],magic[2]);

  // NETCDF
  if (magic[0]==0x43 && magic[1]==0x44 && magic[2]==0x46) {
    if (debug>0) fprintf(stdout,"  NETCDF file\n");
    fileFormat=AMBERNETCDF;
    if (compressType!=NONE) {
      fprintf(stdout,"Error: Compressed NETCDF files are not currently supported.\n");
      return 1;
    }
    return 0;
  }

  // ID by file characteristics; read the first two lines
  IO->Open(filename,"r"); // NOTE: Err Check
  IO->Gets(buffer1,BUFFER_SIZE);
  IO->Gets(buffer2,BUFFER_SIZE);
  IO->Close();

  // Check for terminal CR before newline, indicates DOS file
  i = strlen(buffer1);
  if ( i>1 ) {
    if (buffer1[ i - 2 ] == '\r') {
      if (debug>0) fprintf(stdout,"  [DOS]");
      isDos=1;
    }
  }

  // If both lines have PDB keywords, assume PDB
  if (isPDBkeyword(buffer1) && isPDBkeyword(buffer2)) {
    if (debug>0) fprintf(stdout,"  PDB file\n");
    fileFormat=PDBFILE;
    return 0;
  }

  // If the %VERSION and %FLAG identifiers are present, assume amber parm
  if (strncmp(buffer1,"%VERSION",8)==0 && strncmp(buffer2,"%FLAG",5)==0) {
    if (debug>0) fprintf(stdout,"  AMBER TOPOLOGY file\n");
    fileFormat=AMBERPARM;
    return 0;
  }

  // Amber Restart
  // Check for an integer (I5) followed by 0-2 scientific floats (E15.7)
  if (strlen(buffer2)<=36) {
    //fprintf(stdout,"DEBUG: Checking restart.\n");
    //fprintf(stdout,"DEBUG: buffer2=[%s]\n",buffer2);
    for (i=0; i<5; i++) {
      if (!isspace(buffer2[i]) && !isdigit(buffer2[i])) break;
      //fprintf(stdout,"DEBUG:    %c is a digit/space.\n",buffer2[i]);
    }
    //fprintf(stdout,"DEBUG: i=%i\n");
    //if ( i==5 && strchr(buffer2,'E')!=NULL ) {
    if ( i==5 ) {
      if (debug>0) fprintf(stdout,"  AMBER RESTART file\n");
      fileFormat=AMBERRESTART;
      return 0;
    }
  }

  // If the second line is 81 bytes, assume Amber Traj
  // If second line is 42 bytes, check for Amber Traj with REMD
  // NOTE: Check for digit every 8 char?
  // NOTE: This wont work for trajectories with <3 coords.
  if (strlen(buffer2)==81) {
    //if ( sscanf(buffer2, "%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2,
    //            TrajCoord+3, TrajCoord+4, TrajCoord+5, TrajCoord+6, TrajCoord+7, TrajCoord+8,
    //            TrajCoord+9) == 10 )
    // Make sure we can read at least 3 coords of width 8
    if ( sscanf(buffer2, "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 ) {
      if (debug>0) fprintf(stdout,"  AMBER TRAJECTORY file\n");
      fileFormat=AMBERTRAJ;
      return 0;
    }
  } else if (strlen(buffer2)==42) {
    if ( strncmp(buffer2,"REMD",4)==0 ||
         strncmp(buffer2,"HREMD",5)==0   ) {
      if (debug>0) fprintf(stdout,"  AMBER TRAJECTORY with (H)REMD header.\n");
      fileFormat=AMBERTRAJ;
      return 0;
    }
  }

  // Unidentified file
  fprintf(stdout,"  Warning: %s: UNKNOWN FILE FORMAT.\n",filename);
  return 1; 
}

