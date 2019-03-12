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

#ifdef MPI
int DataSet_Coords_CRD::Sync(size_t total, std::vector<int> const& rank_frames,
                             Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    coords_.resize( total, std::vector<float>( numCrd_+numBoxCrd_ ) );
    int cidx = rank_frames[0]; // Index on master
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int ridx = 0; ridx != rank_frames[rank]; ridx++, cidx++)
        commIn.SendMaster( &(coords_[cidx][0]), numCrd_+numBoxCrd_, rank, MPI_FLOAT );
    }
  } else // Send data to master
    for (unsigned int ridx = 0; ridx != coords_.size(); ++ridx)
      commIn.SendMaster( &(coords_[ridx][0]), numCrd_+numBoxCrd_, commIn.Rank(), MPI_FLOAT );
  return 0;
}
#endif
