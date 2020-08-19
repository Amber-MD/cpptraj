#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

/** Reserve space for coords. */
int DataSet_Coords_CRD::Allocate(SizeArray const& sizeIn) {
  if (!sizeIn.empty())
    coords_.reserve( sizeIn[0] );
  return 0;
}

/** Allocate space in coords_ array. */
int DataSet_Coords_CRD::MemAlloc( SizeArray const& sizeIn ) {
  mprintf("DEBUG: Resize %s to %zu\n", legend(), sizeIn[0]);
  if (!sizeIn.empty()) {
    coords_.resize( sizeIn[0] );
  }
  return 0;
}

/** Copy block from incoming set of same type. */
void DataSet_Coords_CRD::CopyBlock(size_t startIdx, DataSet const* dptrIn, size_t pos, size_t nelts)
{
  DataSet_Coords_CRD const& setIn = static_cast<DataSet_Coords_CRD const&>( *dptrIn );
  // Check that this is compatible
  if (numCrd_ == 0) {
    CoordsSetup( setIn.top_, setIn.cInfo_ );
  } else {
    if (numCrd_ != setIn.numCrd_ || numBoxCrd_ != setIn.numBoxCrd_) {
      mprinterr("Error: Cannot set %s sizes do not match set %s, cannot copy.\n",
                legend(), setIn.legend());
      return;
    }
  }
  CRDarray::iterator begin = coords_.begin() + startIdx;
  std::fill( begin, begin + nelts, std::vector<float>( numCrd_+numBoxCrd_ ) );
  CRDarray::const_iterator ptr = setIn.coords_.begin() + pos;
  std::copy( ptr, ptr + nelts, begin );
}

/** Set up COORDS with given Topology and coordinate info. */
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
