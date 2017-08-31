#include "CoordinateInfo.h"
#include "CpptrajStdio.h"

void CoordinateInfo::PrintCoordInfo(const char* name, const char* parm) const {
  mprintf("DBG: '%s' parm '%s' CoordInfo={ box type %s", name, parm, box_.TypeName());
  if (remdDim_.Ndims() > 0) mprintf(", %i rep dims", remdDim_.Ndims());
  if (hasVel_) mprintf(", velocities");
  if (hasTemp_) mprintf(", temps");
  if (hasTime_) mprintf(", times");
  if (hasFrc_) mprintf(", forces");
  if (ensembleSize_ > 0) mprintf(", ensemble size %i", ensembleSize_);
  mprintf(" }\n");
}

static inline void Append(std::string& meta, std::string const& str) {
  if (meta.empty())
    meta.assign( str );
  else
    meta.append(", " + str);
}

std::string CoordinateInfo::InfoString() const {
  std::string meta;
  if ( HasCrd() )         Append(meta, "coordinates");
  if ( HasVel() )         Append(meta, "velocities");
  if ( HasForce() )       Append(meta, "forces");
  if ( HasTemp() )        Append(meta, "temperature");
  if ( HasTime() )        Append(meta, "time");
  if ( HasReplicaDims() ) Append(meta, "replicaDims");
  if ( HasBox() )         Append(meta, "box");
  return meta;
}

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
