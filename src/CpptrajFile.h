#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
// CpptrajFile
// Class to abstract handling of basic file routines.
// NOTE: Should Format/Type etc stuff be in their own files?
#include "BaseFileIO.h"
// BUFFER_SIZE: Used in Action_NAstruct.cpp FortranFormat.cpp main.cpp 
// CpptrajFile.cpp CpptrajFile.h Traj_AmberCoord.cpp
#define BUFFER_SIZE 1024 
#define OUTPUTFRAMESHIFT 1 // Used for output in DataFile and some TrajFiles 

/* FILE FORMAT:
 * File format specifies how the data in the file is organized. Not used by
 * Ptrajfile itself but by higher-level classes.
 */
enum FileFormat {
  UNKNOWN_FORMAT, PDBFILE, AMBERTRAJ, AMBERNETCDF, AMBERPARM, 
  DATAFILE, AMBERRESTART, AMBERREMD, XMGRACE, CONFLIB, 
  AMBERRESTARTNC, MOL2FILE, GNUPLOT, CHARMMPSF, CHARMMDCD
};

/* FILE TYPE:
 * File type describes how the file is accessed at the lowest level and
 * determines what IO class is associated with the file.
 */
enum FileType {
  UNKNOWN_TYPE, STANDARD, GZIPFILE, BZIP2FILE, ZIPFILE, MPIFILE
};

/* COMPRESS TYPE:
 * Type of compression present if any.
 */
enum CompressType {
  NONE, GZIP, BZIP2, ZIP
};

// ACCESS TYPE
enum AccessType {
  READ, WRITE, APPEND
};

class CpptrajFile {
  protected:
    typedef char enumToken[30];
    static const enumToken FileFormatList[];
    static const enumToken FileTypeList[];
    static const enumToken AccessList[];

    int isOpen;
    int debug;

    void determineType();
    void determineFormat();
    void SetBaseFilename();
    int SetupRead();
    int SetupWrite();
  public:
    AccessType access;
    off_t uncompressed_size;
    off_t file_size;
    BaseFileIO *IO; 
    FileType fileType;
    FileFormat fileFormat;
    char *filename;                  // Passed in filename
    char *basefilename;              // Filename minus any path
    char *Ext;                       // Filename extension
    CompressType compressType;
    int isDos;                       // 1 if CR present, need to count them as newlines

    CpptrajFile();
    ~CpptrajFile();

    FileFormat GetFmtFromArg(char *, FileFormat);
    void SetExtFromFmt(char *, FileFormat);
    char *Type(FileType);
    char *Format(FileFormat);
    int SetupFile(char *, AccessType, int); 
    int SetupFile(char*,AccessType,FileFormat,FileType,int);
    int OpenFile();
    void CloseFile();

    bool IsOpen() { if (isOpen) return true; else return false; }
};
#endif
