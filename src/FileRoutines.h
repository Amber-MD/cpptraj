#ifndef INC_FILEROUTINES_H
#define INC_FILEROUTINES_H
/*! \file FileRoutines.h
    \brief Definitions and routines used for files.
 */
// NOTE: When adding to FileFormat etc be sure to update FileFormatList
// in FileRoutines.c

// FILE FORMAT:
/** File format specifies how the data in the file is organized. Returned 
  * by the file detection routine in CpptrajFile::SetupRead, and by 
  * higher-level classes such as TrajectoryIO, DataFile, and AmberParm.
  */
enum FileFormat {
  UNKNOWN_FORMAT, PDBFILE, AMBERTRAJ, AMBERNETCDF, AMBERPARM,
  DATAFILE, AMBERRESTART, AMBERREMD, XMGRACE, CONFLIB,
  AMBERRESTARTNC, MOL2FILE, GNUPLOT, CHARMMPSF, CHARMMDCD,
  OLDAMBERPARM
};
    
// FILE TYPE:
/** File type describes how the file is accessed at the lowest level and
  * determines what FileIO class is associated with the file.
  */ 
enum FileType {
  UNKNOWN_TYPE, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE
};  
    
// COMPRESS TYPE:
/** Type of compression present if any.  */
enum CompressType {
  NO_COMPRESSION, GZIP, BZIP2, ZIP
};

// ACCESS TYPE
enum AccessType {
  READ, WRITE, APPEND
};

char *File_Format(FileFormat);
char *File_Type(FileType);
char *File_Access(AccessType);
FileFormat GetFmtFromArg(char *, FileFormat);
void SetExtFromFmt(char *, FileFormat);
FileType DetermineType(FileType, char *);
FileFormat DetermineFormat(FileFormat, char *);
#endif
