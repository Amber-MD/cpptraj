#include "DataSet_Coords_FRM.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
DataSet_Coords_FRM::DataSet_Coords_FRM() {}

// -----------------------------------------------
#ifdef MPI
// FIXME Sync
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


