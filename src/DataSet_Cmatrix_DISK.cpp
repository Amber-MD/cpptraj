#include "DataSet_Cmatrix_DISK.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

int DataSet_Cmatrix_DISK::AllocateCmatrix(size_t sizeIn) {
  if (fname_.empty()) {
    mprinterr("Internal Error: Cluster matrix file name not set.\n");
    return 1;
  }
  mprintf("\tEstimated pair-wise matrix disk usage: > %s\n",
          ByteString( ((sizeIn*(sizeIn-1))/2)*sizeof(float), BYTE_DECIMAL).c_str());
  if (file_.OpenCmatrixWrite(fname_, sievedFrames_.MaxFrames(), sizeIn, sievedFrames_.Sieve()))
    return 1;
  // Write actual frames array if necessary
  if (sievedFrames_.Type() != ClusterSieve::NONE) {
    if (file_.WriteFramesArray( sievedFrames_.Frames() ))
      return 1;
  }
  return 0;
}
