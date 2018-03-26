#include "FileIO_MpiShared.h"
#ifdef MPI
// Needed for some older Intel MPI and newer OpenMPI versions
# include <cstdio>

// FileIO_MpiShared::Write()
int FileIO_MpiShared::Write(const void *buffer, size_t num_bytes) {
  if (pfile_.Fwrite_shared(buffer, num_bytes, MPI_CHAR)) return 1;
  // NOTE: Check for errors here.
  return 0;
}
#endif
