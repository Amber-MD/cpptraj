#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"

int DataSet_Coords_CRD::AllocateCoords( size_t sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

double DataSet_Coords_CRD::sizeInMB(size_t nframes, size_t natom, size_t nbox) {
  size_t frame_size_bytes = ((natom * 3UL) + nbox) * sizeof(float);
  double sze = (double)((nframes * frame_size_bytes) + sizeof(CRDarray));
  return (sze / (1024 * 1024));
}

void DataSet_Coords_CRD::Info() const {
  double sze = sizeInMB(coords_.size(), top_.Natom(), numBoxCrd_);
  if (sze < 1.0)
    mprintf(" (<1 MB)");
  else
    mprintf(" (%.2f MB)", sze);
  if (numBoxCrd_ > 0) mprintf(" Box Coords,");
  if (hasVel_)        mprintf(" Velocities,");
  mprintf(" %i atoms", top_.Natom());
}
