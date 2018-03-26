#ifndef INC_FILEIO_MPI_H
#define INC_FILEIO_MPI_H
#ifdef MPI
#include "FileIO.h" 
#include "Parallel.h"
/// MPI file IO
class FileIO_Mpi : public FileIO {
  public:
    FileIO_Mpi() {}
    virtual ~FileIO_Mpi() {}
    int OpenStream(StreamType) { return 1; } 
    int Open(const char *, const char *);    
    int Close();
    int Read(void *, size_t );
    virtual int Write(const void *, size_t);
    int Flush();
    int Seek(off_t);
    int Rewind();  
    off_t Tell();
    int Gets(char *, int );
    int SetSize(long int);
    off_t Size(const char*) { return 0; }
    /// Set communicator
    void SetComm(Parallel::Comm const& commIn) { comm_ = commIn; }
  protected:
    Parallel::File pfile_;
    Parallel::Comm comm_; 
};
#endif
#endif
