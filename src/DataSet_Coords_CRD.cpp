#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString
#include <algorithm> // std::copy

/** CONSTRUCTOR */
DataSet_Coords_CRD::DataSet_Coords_CRD() :
  DataSet_Coords(COORDS),
  framesToReserve_(0)
{}

/** Reserve space for coords. */
int DataSet_Coords_CRD::Allocate(SizeArray const& sizeIn) {
  //rprintf("DEBUG: Calling Allocate: SizeArray[0]= %zu  framesToReserve= %i\n", sizeIn[0], framesToReserve_);
  if (!sizeIn.empty()) {
    framesToReserve_ = (int)sizeIn[0];
    if (frames_.HasComponents())
      frames_.Resize( framesToReserve_ );
  }
  return 0;
}

/** Allocate space in coords_ array. */
int DataSet_Coords_CRD::MemAlloc( SizeArray const& sizeIn ) {
  //rprintf("DEBUG: MemAlloc resize %s to %zu\n", legend(), sizeIn[0]);
  if (!sizeIn.empty()) {
    framesToReserve_ = (int)sizeIn[0];
    frames_.Resize( framesToReserve_ );
  }
  return 0;
}

/** Copy block from incoming set of same type.
  * \param startIdx Position in this set to copy to.
  * \param dptrIn   Set to copy from.
  * \param pos      Position in dptrIn to copy from.
  * \param nelts    Number of elements from dptrIn to copy.
  */
void DataSet_Coords_CRD::CopyBlock(size_t startIdx, DataSet const* dptrIn, size_t pos, size_t nelts)
{
  DataSet_Coords_CRD const& setIn = static_cast<DataSet_Coords_CRD const&>( *dptrIn );
  // Check that this is compatible
  if (!frames_.HasComponents()) {
    CoordsSetup( setIn.top_, setIn.cInfo_ );
  } else {
    if (frames_ != setIn.frames_) {
      mprinterr("Error: Cannot set %s sizes do not match set %s, cannot copy.\n",
                legend(), setIn.legend());
      return;
    }
  }
  CompactFrameArray::iterator begin = frames_.frameBegin(startIdx);
  CompactFrameArray::const_iterator ptr = setIn.frames_.frameBegin(pos);
  std::copy( ptr, ptr + (nelts*frames_.FrameSize()), begin);
}

/** Set up COORDS with given Topology and coordinate info. */
int DataSet_Coords_CRD::CoordsSetup(Topology const& topIn, CoordinateInfo const& cInfoIn) {
  top_ = topIn;
  cInfo_ = cInfoIn;

  if (frames_.SetupFrameArray(cInfo_, topIn.Natom(), framesToReserve_)) {
    mprinterr("Internal Error: Could not set up CompactFrameArray for '%s'\n", legend());
    return 1;
  }

  return 0;
}

/** Convert frame to CompactFrame */
static inline void FrameToArray(CompactFrameArray& frames_, Frame const& fIn) {
  if (frames_.HasComponent(CoordinateInfo::POSITION)) frames_.SetFromDblPtr(fIn.xAddress(), CoordinateInfo::POSITION);
  if (frames_.HasComponent(CoordinateInfo::VELOCITY)) frames_.SetFromDblPtr(fIn.vAddress(), CoordinateInfo::VELOCITY);
  if (frames_.HasComponent(CoordinateInfo::FORCE)) frames_.SetFromDblPtr(fIn.fAddress(), CoordinateInfo::FORCE);
  if (frames_.HasComponent(CoordinateInfo::BOX)) {
    // Prefer storing XYZ ABG if possible; loses less info, especially when converting to float.
    if (fIn.BoxCrd().Is_X_Aligned()) {
      double dtmp[9];
      std::copy(fIn.BoxCrd().XyzPtr(), fIn.BoxCrd().XyzPtr()+6, dtmp);
      std::fill(dtmp+6, dtmp+9, 0);
      frames_.SetFromDblPtr(dtmp, CoordinateInfo::BOX);
    } else {
      frames_.SetFromDblPtr(fIn.BoxCrd().UnitCell().Dptr(), CoordinateInfo::BOX);
    }
  }
  if (frames_.HasComponent(CoordinateInfo::TEMPERATURE)) frames_.SetFromDblVal(fIn.Temperature(), CoordinateInfo::TEMPERATURE);
  if (frames_.HasComponent(CoordinateInfo::PH)) frames_.SetFromDblVal(fIn.pH(), CoordinateInfo::PH);
  if (frames_.HasComponent(CoordinateInfo::REDOX)) frames_.SetFromDblVal(fIn.RedOx(), CoordinateInfo::REDOX);
  if (frames_.HasComponent(CoordinateInfo::TIME)) frames_.SetFromDblVal(fIn.Time(), CoordinateInfo::TIME);
  if (frames_.HasComponent(CoordinateInfo::STEP)) frames_.SetFromIntVal(fIn.Step(), CoordinateInfo::STEP);
  if (frames_.HasComponent(CoordinateInfo::REMD_INDICES)) frames_.SetFromIntPtr(fIn.iAddress(), CoordinateInfo::REMD_INDICES);
  if (frames_.HasComponent(CoordinateInfo::REPIDX)) frames_.SetFromIntVal(fIn.RepIdx(), CoordinateInfo::REPIDX);
  if (frames_.HasComponent(CoordinateInfo::CRDIDX)) frames_.SetFromIntVal(fIn.CrdIdx(), CoordinateInfo::CRDIDX);
}

/** Add frame to array. */
void DataSet_Coords_CRD::AddFrame(Frame const& fIn) {
  frames_.NextAndAllocate();

  FrameToArray(frames_, fIn);
}

/** Copy frame to specified position in array. */
void DataSet_Coords_CRD::SetCRD(int idx, Frame const& fIn) {
  frames_.SeekAndAllocate(idx);

  FrameToArray(frames_, fIn);
}

/** Non-coord-related array to Frame */
static inline void ArrayToFrame(int idx, Frame& fOut, CompactFrameArray const& frames_) {
  // TODO should not have to check incoming box
  if (fOut.BoxCrd().HasBox() && frames_.HasComponent(CoordinateInfo::BOX)) {
    double dtmp[9];
    frames_.GetToDblPtr(dtmp, idx, CoordinateInfo::BOX);
    // Box may be stored as XYZ ABG or unit cell vectors. If last three elements are all zero
    // it is the former, otherwise the latter.
    if (dtmp[6] != 0 || dtmp[7] != 0 || dtmp[8] != 0)
      fOut.ModifyBox().AssignFromUcell( dtmp );
    else
      fOut.ModifyBox().AssignFromXyzAbg( dtmp );
  }
  if (frames_.HasComponent(CoordinateInfo::TEMPERATURE)) fOut.SetTemperature(frames_.GetVal(idx, CoordinateInfo::TEMPERATURE));
  if (frames_.HasComponent(CoordinateInfo::PH)) fOut.Set_pH(frames_.GetVal(idx, CoordinateInfo::PH));
  if (frames_.HasComponent(CoordinateInfo::REDOX)) fOut.SetRedOx(frames_.GetVal(idx, CoordinateInfo::REDOX));
  if (frames_.HasComponent(CoordinateInfo::TIME)) fOut.SetTime(frames_.GetVal(idx, CoordinateInfo::TIME));
  if (frames_.HasComponent(CoordinateInfo::STEP)) fOut.SetStep(frames_.GetVal(idx, CoordinateInfo::STEP));
  if (frames_.HasComponent(CoordinateInfo::REMD_INDICES) && !fOut.RemdIndices().empty()) frames_.GetToIntPtr(fOut.iAddress(), idx, CoordinateInfo::REMD_INDICES);
  if (frames_.HasComponent(CoordinateInfo::REPIDX)) fOut.SetRepIdx(frames_.GetVal(idx, CoordinateInfo::REPIDX));
  if (frames_.HasComponent(CoordinateInfo::CRDIDX)) fOut.SetCrdIdx(frames_.GetVal(idx, CoordinateInfo::CRDIDX));
}

/** Get a frame from specified position in array. */
void DataSet_Coords_CRD::GetFrame(int idx, Frame& fOut) {
  if (frames_.HasComponent(CoordinateInfo::POSITION)) frames_.GetToDblPtr(fOut.xAddress(), idx, CoordinateInfo::POSITION);
  if (frames_.HasComponent(CoordinateInfo::VELOCITY) && fOut.HasVelocity()) frames_.GetToDblPtr(fOut.vAddress(), idx, CoordinateInfo::VELOCITY);
  if (frames_.HasComponent(CoordinateInfo::FORCE) && fOut.HasForce()) frames_.GetToDblPtr(fOut.fAddress(), idx, CoordinateInfo::FORCE);
  ArrayToFrame(idx, fOut, frames_);
}

/** Get selected atoms from a frame from specified position in array. */
void DataSet_Coords_CRD::GetFrame(int idx, Frame& fOut, AtomMask const& mask) {
  if (frames_.HasComponent(CoordinateInfo::POSITION)) frames_.GetToMaskDblPtr(fOut.xAddress(), mask.Selected(), idx, CoordinateInfo::POSITION);
  if (frames_.HasComponent(CoordinateInfo::VELOCITY) && fOut.HasVelocity()) frames_.GetToMaskDblPtr(fOut.vAddress(), mask.Selected(), idx, CoordinateInfo::VELOCITY);
  if (frames_.HasComponent(CoordinateInfo::FORCE) && fOut.HasForce()) frames_.GetToMaskDblPtr(fOut.fAddress(), mask.Selected(), idx, CoordinateInfo::FORCE);
  ArrayToFrame(idx, fOut, frames_);
}

/*
size_t DataSet_Coords_CRD::sizeInBytes(size_t nframes, size_t natom, size_t nbox) {
  size_t frame_size_bytes = ((natom * 3UL) + nbox) * sizeof(float);
  return ((nframes * frame_size_bytes) + sizeof(CRDarray));
}*/

#ifdef MPI
int DataSet_Coords_CRD::Sync(size_t total, std::vector<int> const& rank_frames,
                             Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  if (commIn.Master()) {
    // Resize for total number of frames.
    frames_.Resize( total );
    int cidx = rank_frames[0]; // Frame Index on master
    // Receive data from each rank
    for (int rank = 1; rank < commIn.Size(); rank++) {
      //rprintf("DEBUG: MASTER RECEIVE from rank %i cidx=%i rank_frames=%i FrameSize=%u\n",
      //        rank, cidx, rank_frames[rank], frames_.FrameSize());
      commIn.SendMaster( frames_.FramePtr(cidx), rank_frames[rank]*frames_.FrameSize(),
                         rank, MPI_FLOAT );
      cidx += rank_frames[rank];
    }
  } else {// Send data to master
    //rprintf("DEBUG: RANK %i SEND TO MASTER Size=%zu FrameSize=%u\n",
    //        commIn.Rank(), Size(), frames_.FrameSize());
    commIn.SendMaster( frames_.FramePtr(0), Size()*frames_.FrameSize(), 
                       commIn.Rank(), MPI_FLOAT );
  }
  return 0;
}
#endif
