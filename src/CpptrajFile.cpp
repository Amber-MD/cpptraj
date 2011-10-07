// CpptrajFile 
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <sys/stat.h>
// NOTE: It seems some PGI compilers do not function correctly when glob.h
//       is included and large file flags are set. Just disable globbing
//       for PGI.
#ifndef __PGI
#  include <glob.h> // For tilde expansion
#endif
#include "CpptrajFile.h"
#include "NetcdfRoutines.h"
#include "PDBfileRoutines.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// File Types
#include "StdFile.h"
#ifdef HASGZ
#  include "GzipFile.h"
#endif
#ifdef MPI
#  include "MpiFile.h"
#endif
#ifdef HASBZ2
#  include "Bzip2File.h"
#endif

//typedef char enumToken[30];
const CpptrajFile::enumToken CpptrajFile::FileFormatList[15] = {
  "UNKNOWN_FORMAT", "PDBFILE", "AMBERTRAJ", "AMBERNETCDF", "AMBERPARM", 
  "DATAFILE", "AMBERRESTART", "AMBERREMD", "XMGRACE", "CONFLIB", 
  "AMBERRESTARTNC", "MOL2FILE", "GNUPLOT", "CHARMMPSF", "CHARMMDCD"
};
const CpptrajFile::enumToken CpptrajFile::FileTypeList[6] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};
const CpptrajFile::enumToken CpptrajFile::AccessList[3] = {
  "R", "W", "A"
};

// CONSTRUCTOR
CpptrajFile::CpptrajFile() {
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
char *CpptrajFile::Format(FileFormat FFin) {
  return (char*)FileFormatList[FFin];
}

// Return enumerated type for FileType
char *CpptrajFile::Type(FileType FTin) {
  return (char*)FileTypeList[FTin];
}

// DESTRUCTOR
CpptrajFile::~CpptrajFile() {
   //fprintf(stderr,"CPPTRAJFILE DESTRUCTOR\n");
   CloseFile();
   delete IO;
   if (filename!=NULL) free(filename);
   if (basefilename!=NULL) free(basefilename);
   if (Ext!=NULL) free(Ext);
}

/* CpptrajFile::GetFmtFromArg()
 * Return file format given a file format keyword. Default to def. 
 * NOTE: def should probably not be allowed to be UNKNOWN_FORMAT,
 * but this is currently not explicitly checked.
 */
FileFormat CpptrajFile::GetFmtFromArg(char *argIn, FileFormat def) {
  FileFormat writeFormat = def;
  if (argIn==NULL) return writeFormat;
  if      ( strcmp(argIn,"pdb")==0      ) writeFormat=PDBFILE;
  else if ( strcmp(argIn,"data")==0     ) writeFormat=DATAFILE;
  else if ( strcmp(argIn,"netcdf")==0   ) writeFormat=AMBERNETCDF;
  else if ( strcmp(argIn,"restart")==0  ) writeFormat=AMBERRESTART;
  else if ( strcmp(argIn,"ncrestart")==0) writeFormat=AMBERRESTARTNC;
  else if ( strcmp(argIn,"restartnc")==0) writeFormat=AMBERRESTARTNC;
  else if ( strcmp(argIn,"mol2")==0     ) writeFormat=MOL2FILE;
  return writeFormat;
}

/* CpptrajFile::SetExtFromFmt()
 * Set buffer with a filename extension corresponding to the given file 
 * format. Currently the longest extension requires buffer size of 7.
 */
void CpptrajFile::SetExtFromFmt(char *buffer, FileFormat fmtIn) {
  if (buffer==NULL) return;
  switch (fmtIn) {
    case PDBFILE       : strcpy(buffer,".pdb"  ); break;
    case AMBERTRAJ     : strcpy(buffer,".crd"  ); break;
    case AMBERNETCDF   : strcpy(buffer,".nc"   ); break;
    case AMBERPARM     : strcpy(buffer,".parm7"); break;
    case DATAFILE      : strcpy(buffer,".dat"  ); break;
    case AMBERRESTART  : strcpy(buffer,".rst7" ); break;
    case XMGRACE       : strcpy(buffer,".agr"  ); break;
    case AMBERRESTARTNC: strcpy(buffer,".ncrst"); break;
    case MOL2FILE      : strcpy(buffer,".mol2" ); break;
    case GNUPLOT       : strcpy(buffer,".gnu"  ); break;
    case CHARMMPSF     : strcpy(buffer,".psf"  ); break;
    case CHARMMDCD     : strcpy(buffer,".dcd"  ); break;
    default:
      strcpy(buffer,"");
  }
}

/* CpptrajFile::OpenFile()
 * Open the file. If already open, reopen.
 */
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
      
  isOpen=1;
  return 0;
}

/* CpptrajFile::CloseFile()
 * Close the file.
 */
void CpptrajFile::CloseFile() {
  if (isOpen) {
    IO->Close();
    if (debug>0) rprintf("Closed %s.\n",filename);
    isOpen=0;
  }
}

/* CpptrajFile::determineType()
 * If the file type is unknown attempt to determine it from filename extension.
 * Default to standard.
 */
void CpptrajFile::determineType() {

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

/* CpptrajFile::determineFormat()
 * If the file format is unknown attempt to determine from filename extension.
 * Default to datafile.
 */
void CpptrajFile::determineFormat() {

  if (fileFormat!=UNKNOWN_FORMAT) return;

  // No extension
  if (Ext==NULL) return;

  if      ( strcmp(Ext,".dat")==0 ) fileFormat=DATAFILE;
  else if ( strcmp(Ext,".agr")==0 ) fileFormat=XMGRACE;
  else if ( strcmp(Ext,".gnu")==0 ) fileFormat=GNUPLOT;
  else fileFormat=DATAFILE;
}  

/* CpptrajFile::SetBaseFilename()
 * Strip leading path from input filename. Use strtok routine to separate 
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

/* CpptrajFile::SetupFile()
 * Set up the given file for the specified access, assume DATAFILE format
 * (no autodetection of format). Implicitly requests autodetection of the 
 * file type.
 */
int CpptrajFile::SetupFile(char *filenameIn, AccessType accessIn, int debugIn) {
  return SetupFile(filenameIn,accessIn,DATAFILE,UNKNOWN_TYPE,debugIn);
}

/* CpptrajFile::SetupFile()
 * Sets the file name, access type (R/W/A), file type and file format (for 
 * WRITE), and debug level. If called with READ or append file type and format
 * will be determined by SetupRead. If called with WRITE the given type and 
 * format will be used; if called with UNKNOWN type and format they will
 * try to be determined by the given file extension. 
 * Can be called with NULL for write, this will write to stdout.
 */
int CpptrajFile::SetupFile(char *filenameIn, AccessType accessIn, 
                         FileFormat fileFormatIn, FileType fileTypeIn, 
                         int debugIn) {
#ifndef __PGI
  glob_t globbuf;
#endif
  debug=debugIn;
  if (debug>0) {
    mprintf("CpptrajFile: [%s] FMT %s, TYPE %s, ACC %s\n",filenameIn,
           FileFormatList[fileFormatIn],FileTypeList[fileTypeIn],
           AccessList[accessIn]);
  }

  // Store filename
  if (filenameIn!=NULL) {
    filename=(char*) malloc( (strlen(filenameIn)+1) * sizeof(char) );
    strcpy(filename,filenameIn);
  }
#ifndef __PGI
  // On read or append do tilde expansion and store new filename
  if (accessIn!=WRITE) {
    // If no filename this is an error.
    if (filename==NULL) {
      mprintf("Error: CpptrajFile::SetupFile: NULL filename specified on READ or APPEND\n");
      return 1;
    }
    globbuf.gl_offs = 1;
    if ( glob(filename, GLOB_TILDE, NULL, &globbuf)!=0 ) return 1; 
    if (debug>1) mprintf("\tGLOB(0): [%s]\n",globbuf.gl_pathv[0]);
    filename=(char*) realloc( filename, (strlen(globbuf.gl_pathv[0])+1) * sizeof(char));
    strcpy(filename, globbuf.gl_pathv[0]);
    globfree(&globbuf);
  }
#endif
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
      fileType=fileTypeIn;
      fileFormat=fileFormatIn;
      this->determineType();
      this->determineFormat();
      if ( SetupWrite() ) return 1; 
      break;
    default: return 1;
  }
  if (debug>0)
    rprintf("\t[%s] is format %s and type %s with access %s\n",filename,
           FileFormatList[fileFormat],FileTypeList[fileType],AccessList[access]);
  return 0;
}

/* CpptrajFile::SetupWrite()
 * Set up file with specified type for writing. fileFormat is set by SetupFile.
 * NOTE: Determine if compression is requested, either by arg or from name
 */
int CpptrajFile::SetupWrite() {

  if (debug>1) mprintf("DEBUG: Setting up WRITE for file %s\n",filename);
  // Eventually allow other file types
  //fileType=STANDARD;
  switch (fileType) {
    case GZIPFILE  : 
#ifdef HASGZ
      IO = new GzipFile(); 
#else
      mprintf("Error: SetupWrite(%s):\n",filename);
      mprintf("       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE :
#ifdef HASBZ2 
      IO = new Bzip2File();
#else
      mprintf("Error: SetupWrite(%s):\n",filename);
      mprintf("       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
      return 1;
#endif
    break;
    //case ZIPFILE   : IO = new ZipFile(); break;
    case STANDARD  : IO = new StdFile();  break;
    case MPIFILE   : 
#ifdef MPI
      IO = new MpiFile();
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

/* CpptrajFile::SetupRead() 
 * Open the file specified by filename for READ or APPEND access. Attempt to 
 * identify the file type and format.
 */
int CpptrajFile::SetupRead() {
  unsigned char magic[3];
  char buffer1[BUFFER_SIZE], buffer2[BUFFER_SIZE];
  char *CheckConventions; // Only used to check if netcdf is traj or restart
  float TrajCoord[10];
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
  IO = new StdFile();

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
  if (access==APPEND && compressType!=NONE) {
    mprintf("Error: Appending to compressed files is not supported.\n");
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
      mprintf("Error: SetupRead(%s):\n",filename);
      mprintf("       Compiled without Gzip support. Recompile with -DHASGZ\n");
      return 1;
#endif
      break;
    case BZIP2FILE : 
#ifdef HASBZ2
      IO = new Bzip2File(); 
#else
      mprintf("Error: SetupRead(%s):\n",filename);
      mprintf("       Compiled without Bzip2 support. Recompile with -DHASBZ2\n");
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

  // If file format already determined, exit.
  if (fileFormat!=UNKNOWN_FORMAT) {
    if (debug>0) mprintf("CpptrajFile: %s: Not detecting file format.\n",filename);
    return 0;
  }

  // ========== Determine format if UNKNOWN ==========
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
    if (compressType!=NONE) {
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

  // If the second line is 81 bytes, assume Amber Traj
  // If second line is 42 bytes, check for Amber Traj with REMD
  // NOTE: Check for digit every 8 char?
  // NOTE: This wont work for trajectories with <3 coords.
  if ((int)strlen(buffer2)==81+isDos) {
    //if ( sscanf(buffer2, "%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2,
    //            TrajCoord+3, TrajCoord+4, TrajCoord+5, TrajCoord+6, TrajCoord+7, TrajCoord+8,
    //            TrajCoord+9) == 10 )
    // Make sure we can read at least 3 coords of width 8
    if ( sscanf(buffer2, "%8f%8f%8f", TrajCoord, TrajCoord+1, TrajCoord+2) == 3 ) {
      if (debug>0) mprintf("  AMBER TRAJECTORY file\n");
      fileFormat=AMBERTRAJ;
      return 0;
    }
  } else if (strlen(buffer2)==42) {
    if ( strncmp(buffer2,"REMD",4)==0 ||
         strncmp(buffer2,"HREMD",5)==0   ) {
      if (debug>0) mprintf("  AMBER TRAJECTORY with (H)REMD header.\n");
      fileFormat=AMBERTRAJ;
      return 0;
    }
  }

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

