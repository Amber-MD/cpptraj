#ifndef INC_CPPTRAJFILE_H
#define INC_CPPTRAJFILE_H
#include "FileIO.h"
#include "FileRoutines.h"
#ifdef USE_CHARBUFFER
#include "CharBuffer.h"
#endif
/// Used in Action_NAstruct.cpp FortranFormat.cpp main.cpp 
/// CpptrajFile.cpp CpptrajFile.h Traj_AmberCoord.cpp
#define BUFFER_SIZE 1024 
/* Compiler Defines:
 * - USE_CHARBUFFER: Use CharBuffer to buffer an entire file.
 */ 
// Class: CpptrajFile
/// Class to abstract handling of basic file routines.
class CpptrajFile {
    bool isOpen;
    int debug;
#   ifdef USE_CHARBUFFER
    CharBuffer c_buffer;
#   endif
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
    char *filename;                  ///< Passed in filename
    char *basefilename;              ///< Filename minus any path
    char *Ext;                       ///< Filename extension
    CompressType compressType;
    int isDos;                       ///< 1 if CR present, need to count them as newlines

    CpptrajFile();
    ~CpptrajFile();

    int SetupFile(char *, AccessType, int); 
    int SetupFile(char*,AccessType,FileFormat,FileType,int);
    int OpenFile();
    void CloseFile();
#   ifdef USE_CHARBUFFER
    int OpenFileBuffered();
    //int ReadBuffered();
    int Gets(char*, int);
    void Rewind();
    int Read(void*,size_t);
#   endif
    bool IsOpen() { if (isOpen) return true; else return false; }
};
#endif
