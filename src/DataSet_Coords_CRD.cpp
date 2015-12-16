#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

int DataSet_Coords_CRD::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty())
    coords_.reserve( sizeIn[0] );
  return 0;
}

int DataSet_Coords_CRD::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;
  numCrd_ = top_.Natom() * 3;
  if (cInfo_.TrajBox().HasBox())
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

size_t DataSet_Coords_CRD::sizeInBytes(size_t nframes, size_t natom, size_t nbox) {
  size_t frame_size_bytes = ((natom * 3UL) + nbox) * sizeof(float);
  return ((nframes * frame_size_bytes) + sizeof(CRDarray));
}

void DataSet_Coords_CRD::Info() const {
  mprintf(" (%s)",
          ByteString(sizeInBytes(coords_.size(), top_.Natom(), numBoxCrd_), BYTE_DECIMAL).c_str());
  CommonInfo();
}
