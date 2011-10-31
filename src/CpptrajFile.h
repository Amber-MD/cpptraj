#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
/// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
#include "FileIO.h"
#include "FileRoutines.h"
// BUFFER_SIZE: Used in Action_NAstruct.cpp FortranFormat.cpp main.cpp 
// CpptrajFile.cpp CpptrajFile.h Traj_AmberCoord.cpp
#define BUFFER_SIZE 1024 
#define OUTPUTFRAMESHIFT 1 // Used for output in DataFile and some TrajFiles 
class CpptrajFile {
    int isOpen;
    int debug;

    void SetBaseFilename();
    int SetupRead();
    int SetupWrite();
  public:
    AccessType access;
    off_t uncompressed_size;
    off_t file_size;
    FileIO *IO; 
    FileType fileType;
    FileFormat fileFormat;
    char *filename;                  // Passed in filename
    char *basefilename;              // Filename minus any path
    char *Ext;                       // Filename extension
    CompressType compressType;
    int isDos;                       // 1 if CR present, need to count them as newlines

    CpptrajFile();
    ~CpptrajFile();

    int SetupFile(char *, AccessType, int); 
    int SetupFile(char*,AccessType,FileFormat,FileType,int);
    int OpenFile();
    void CloseFile();

    bool IsOpen() { if (isOpen) return true; else return false; }
};
#endif
