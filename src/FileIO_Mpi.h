#ifndef INC_FILEIO_MPI_H
#define INC_FILEIO_MPI_H
#include "FileIO.h" 
#include "MpiRoutines.h"
// Class: FileIO_Mpi
/// MPI file IO, wrappers for the MPI routines in MpiRoutines.h
class FileIO_Mpi : public FileIO {
    parallelType pfile; 
  public:
    FileIO_Mpi(); 
    ~FileIO_Mpi(); 

    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t, size_t );
    int Write(void *, size_t, size_t);  
    int Seek(off_t);
    int Rewind();  
    //long int Tell();  // NOTE: Tell may be unnecessary if only for size reporting.
    int Gets(char *, int );
    int SetSize(long int);
};
#endif
