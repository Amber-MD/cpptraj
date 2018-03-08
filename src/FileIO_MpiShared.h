#ifndef INC_FILEIO_MPISHARED_H
#define INC_FILEIO_MPISHARED_H
#ifdef MPI
#include "FileIO_Mpi.h"
class FileIO_MpiShared : public FileIO_Mpi {
  public:
    FileIO_MpiShared() {}
    /// This is a shared write
    int Write(const void*, size_t);
};
#endif
#endif
