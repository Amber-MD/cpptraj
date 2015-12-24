#include "CoordinateInfo.h"
#ifdef MPI
int CoordinateInfo::SyncCoordInfo(Parallel::Comm const& commIn) {
  // ensSize, hasvel, hastemp, hastime, hasfrc, NrepDims, Dim1, ..., DimN, 
  int* iArray;
  int iSize;
  if (commIn.Master()) {
    iSize = remdDim_.Ndims() + 6;
    commIn.MasterBcast( &iSize, 1, MPI_INT );
    iArray = new int[ iSize ];
    iArray[0] = ensembleSize_;
    iArray[1] = (int)hasVel_;
    iArray[2] = (int)hasTemp_;
    iArray[3] = (int)hasTime_;
    iArray[4] = (int)hasFrc_;
    iArray[5] = remdDim_.Ndims();
    unsigned int ii = 6;
    for (int ir = 0; ir != remdDim_.Ndims(); ir++, ii++)
      iArray[ii] = remdDim_[ir];
    commIn.MasterBcast( iArray, iSize, MPI_INT );
  } else {
    commIn.MasterBcast( &iSize, 1, MPI_INT );
    iArray = new int[ iSize ];
    commIn.MasterBcast( iArray, iSize, MPI_INT );
    ensembleSize_ = iArray[0];
    hasVel_       = (bool)iArray[1];
    hasTemp_      = (bool)iArray[2];
    hasTime_      = (bool)iArray[3];
    hasFrc_       = (bool)iArray[4];
    remdDim_.clear();
    unsigned int ii = 6;
    for (int ir = 0; ir != iArray[5]; ir++, ii++)
      remdDim_.AddRemdDimension( iArray[ii] );
  }
  delete[] iArray;
  box_.SyncBox( commIn );
  return 0;
}
#endif
