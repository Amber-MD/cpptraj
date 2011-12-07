#include <cstring>
#include "FileRoutines.h"
// FileRoutines

// When adding to FileFormatList etc be sure to increment the first # of
// the array and update FileFormat in FileRoutines.h 
const char FileFormatList[16][30] = {
  "UNKNOWN_FORMAT", "PDBFILE", "AMBERTRAJ", "AMBERNETCDF", "AMBERPARM",
  "DATAFILE", "AMBERRESTART", "AMBERREMD", "XMGRACE", "CONFLIB",
  "AMBERRESTARTNC", "MOL2FILE", "GNUPLOT", "CHARMMPSF", "CHARMMDCD",
  "OLDAMBERPARM"
};
const char FileTypeList[6][30] = {
  "UNKNOWN_TYPE", "STANDARD", "GZIPFILE", "BZIP2FILE", "ZIPFILE", "MPIFILE"
};
const char AccessList[3][2] = {
  "R", "W", "A"
};

/// Return enumerated type for FileFormat
char *File_Format(FileFormat FFin) {
  return (char*)FileFormatList[FFin];
}

/// Return enumerated type for FileType
char *File_Type(FileType FTin) {
  return (char*)FileTypeList[FTin];
}

/// Return enumerated type for AccessType
char *File_Access(AccessType FAin) {
  return (char*)AccessList[FAin];
}

// GetFmtFromArg()
/** Return file format given a file format keyword. Default to def. 
  * NOTE: def should probably not be allowed to be UNKNOWN_FORMAT,
  * but this is currently not explicitly checked.
  */
FileFormat GetFmtFromArg(char *argIn, FileFormat def) {
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

// SetExtFromFmt()
/** Set buffer with a filename extension corresponding to the given file 
  * format. Currently the longest extension requires buffer size of 7.
  */
void SetExtFromFmt(char *buffer, FileFormat fmtIn) {
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

// DetermineType()
/** If the file type is unknown attempt to determine it from filename extension.
  * Default to standard.
  */
FileType DetermineType(FileType fileType, char *Ext) {
  // If type is not unknown return input type
  if (fileType!=UNKNOWN_TYPE) return fileType;
  // No extension - return standard
  if (Ext==NULL) return STANDARD;
  // Check extension  
  if      ( strcmp(Ext,".gz")==0  ) return GZIPFILE;
  else if ( strcmp(Ext,".bz2")==0 ) return BZIP2FILE;
  // Otherwise return standard 
  return STANDARD;
}

// DetermineFormat()
/** If the file format is unknown attempt to determine from filename extension.
  * Default to datafile.
  */
FileFormat DetermineFormat(FileFormat fileFormat, char *Ext) {
  // If format is not unknown return input format
  if (fileFormat!=UNKNOWN_FORMAT) return fileFormat;
  // No extension - return datafile 
  if (Ext==NULL) return DATAFILE;
  // Check extension
  if      ( strcmp(Ext,".dat")==0 ) return DATAFILE;
  else if ( strcmp(Ext,".agr")==0 ) return XMGRACE;
  else if ( strcmp(Ext,".gnu")==0 ) return GNUPLOT;
  // Otherwise return datafile
  return DATAFILE;
}  

