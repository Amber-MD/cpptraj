#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"

int DataSet_Coords_CRD::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty())
    coords_.reserve( sizeIn[0] );
  return 0;
}

int DataSet_Coords_CRD::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;
  numCrd_ = top_.Natom() * 3;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
  // FIXME: The COORDS DataSet cannot store things like rep dims, times, or
  //        temperatures. Remove these from the CoordinateInfo and warn.
  if (cInfo_.ReplicaDimensions().Ndims() > 0) {
    mprintf("Warning: COORDS data sets do not store replica dimensions.\n");
    cInfo_.SetReplicaDims( ReplicaDimArray() );
  }
  if (cInfo_.HasTemp()) {
    mprintf("Warning: COORDS data sets do not store temperatures.\n");
    cInfo_.SetTemperature( false );
  }
  if (cInfo_.HasTime()) {
    mprintf("Warning: COORDS data sets do not store times.\n");
    cInfo_.SetTime( false );
  }
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
  CommonInfo();
}
