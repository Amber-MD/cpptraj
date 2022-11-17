#include "DataSet_Coords_FRM.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
DataSet_Coords_FRM::DataSet_Coords_FRM() :
  DataSet_Coords(FRAMES)
{}

// -----------------------------------------------
#ifdef MPI
/** Sync all frames to the master process. */
int DataSet_Coords_FRM::Sync(size_t total, std::vector<int> const& rank_frames,
                             Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    frames_.resize( total, AllocateFrame() );
    int cidx = rank_frames[0]; // Frame Index on master
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int idx = 0; idx < rank_frames[rank]; idx++) {
        //rprintf("DEBUG: MASTER RECEIVE from rank %i cidx=%i rank_frames=%i FrameSize=%u\n",
        //        rank, cidx, rank_frames[rank], frames_.FrameSize());
        frames_[cidx++].RecvFrame( rank, commIn );
      }
    }
  } else {// Send data to master
    for (unsigned int idx = 0; idx < Size(); idx++) {
      //rprintf("DEBUG: RANK %i SEND TO MASTER Size=%zu FrameSize=%u\n",
      //        commIn.Rank(), Size(), frames_.FrameSize());
      frames_[idx].SendFrame(0, commIn);
    }
  }
  return 0;
}
#endif

/** Add a single element. */
void DataSet_Coords_FRM::Add(size_t idx, const void* ptrIn) {
  Frame const* frmPtr = static_cast<Frame const*>( ptrIn );
  if (frmPtr->Natom() != top_.Natom()) {
    mprintf("Warning: DataSet_Coords_FRM::Add: Incoming frame has %i atoms but topology '%s' has %i\n",
            frmPtr->Natom(), top_.c_str(), top_.Natom());
  }
  if (idx > frames_.size())
    frames_.resize( idx );
  // Insert at end
  frames_.push_back( *frmPtr );
}

/** Reserve space in frames_ array. */
int DataSet_Coords_FRM::Allocate(SizeArray const& sizeIn) {
  //rprintf("DEBUG: Calling Allocate: SizeArray[0]= %zu  framesToReserve= %i\n", sizeIn[0], framesToReserve_);
  if (!sizeIn.empty()) {
    /*framesToReserve_ = (int)sizeIn[0];
    if (frames_.HasComponents())
      frames_.Resize( framesToReserve_ );*/
    frames_.reserve( sizeIn[0] );
  }
  return 0;
}

/** Allocate space in frames_ array. */
int DataSet_Coords_FRM::MemAlloc( SizeArray const& sizeIn ) {
  //rprintf("DEBUG: MemAlloc resize %s to %zu\n", legend(), sizeIn[0]);
  if (!sizeIn.empty()) {
    //framesToReserve_ = (int)sizeIn[0];
    //frames_.Resize( framesToReserve_ );
    frames_.resize( sizeIn[0] );
  }
  return 0;
}

/** \return Sum of bytes used by all frames in array. */
size_t DataSet_Coords_FRM::MemUsageInBytes() const {
  size_t total = 0;
  for (FrmArrayType::const_iterator it = frames_.begin(); it != frames_.end(); ++it)
    total += it->DataSize();
  return total;
}

/** Copy block from incoming set of same type.
  * \param startIdx Position in this set to copy to.
  * \param dptrIn   Set to copy from.
  * \param pos      Position in dptrIn to copy from.
  * \param nelts    Number of elements from dptrIn to copy.
  */
void DataSet_Coords_FRM::CopyBlock(size_t startIdx, DataSet const* dptrIn, size_t pos, size_t nelts)
{
  DataSet_Coords_FRM const& setIn = static_cast<DataSet_Coords_FRM const&>( *dptrIn );
  size_t tgtidx = startIdx;
  size_t srcend = pos + nelts;
  for (size_t srcidx = pos; srcidx != srcend; srcidx++, tgtidx++)
    frames_[tgtidx] = setIn.frames_[srcidx];
}

// -----------------------------------------------
/** Set up COORDS with given Topology and coordinate info. */
int DataSet_Coords_FRM::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;

  //if (frames_.SetupFrameArray(cInfo_, topIn.Natom(), framesToReserve_)) {
  //  mprinterr("Internal Error: Could not set up CompactFrameArray for '%s'\n", legend());
  //  return 1;
  //}

  return 0;
}

/** Add frame to array. */
void DataSet_Coords_FRM::AddFrame(Frame const& fIn) {
  if (fIn.Natom() != top_.Natom()) {
    mprintf("Warning: DataSet_Coords_FRM::AddFrame: Incoming frame has %i atoms but topology '%s' has %i\n",
            fIn.Natom(), top_.c_str(), top_.Natom());
  }
  frames_.push_back( fIn );
}

/** Copy frame to specified position in array. */
void DataSet_Coords_FRM::SetCRD(int idx, Frame const& fIn) {
  if (fIn.Natom() != top_.Natom()) {
    mprintf("Warning: DataSet_Coords_FRM::SetCrd: Incoming frame has %i atoms but topology '%s' has %i\n",
            fIn.Natom(), top_.c_str(), top_.Natom());
  }
  if (idx > (int)frames_.size()) // TODO fix this sign issue
    frames_.resize( idx );
  frames_[idx] = fIn; 
}

/** Get a frame from specified position in array. */
void DataSet_Coords_FRM::GetFrame(int idx, Frame& fOut) {
  fOut.SetFrame( frames_[idx] );
}

/** Get selected atoms from a frame from specified position in array. */
void DataSet_Coords_FRM::GetFrame(int idx, Frame& fOut, AtomMask const& mask) {
  fOut.SetFrame( frames_[idx], mask );
}
