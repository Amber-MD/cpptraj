#ifndef INC_PTRAJFILE_H
#define INC_PTRAJFILE_H
// PtrajFile
// Class to abstract handling of basic file routines.
#include "BaseFileIO.h"

#define BUFFER_SIZE 1024 // Used in PtrajState, AmberTraj, and PtrajFile
#define OUTPUTFRAMESHIFT 1 // Used for output in DataFile and some TrajFiles 

/* FILE FORMAT:
 * File format specifies how the data in the file is organized. Not used by
 * Ptrajfile itself but by higher-level classes.
 */
enum FileFormat {
  UNKNOWN_FORMAT, PDBFILE, AMBERTRAJ, AMBERNETCDF, AMBERPARM, 
  DATAFILE, AMBERRESTART, AMBERREMD, XMGRACE, CONFLIB, AMBERRESTARTNC,
  MOL2FILE
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

class PtrajFile {
  protected:
    typedef char enumToken[30];
    static const enumToken FileFormatList[];
    static const enumToken FileTypeList[];

    int isOpen;
    int debug;

    void determineType();
    void determineFormat();
    void SetBaseFilename();
    int SetupRead();
    int SetupWrite();
  public:
    AccessType access;
    long long int uncompressed_size;
    off_t file_size;
    BaseFileIO *IO; 
    FileType fileType;
    FileFormat fileFormat;
    char *filename;                  // Passed in filename
    char *basefilename;              // Filename minus any path
    char *Ext;                       // Filename extension
    CompressType compressType;
    int isDos;                       // 1 if CR present, need to count them as newlines

    PtrajFile();
    ~PtrajFile();

    char *Type();
    char *Format();    
    int SetupFile(char*,AccessType,FileFormat,FileType,int);
    int OpenFile();
    void CloseFile();
};
#endif
