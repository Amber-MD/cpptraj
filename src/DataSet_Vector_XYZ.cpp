#include "DataSet_Vector_XYZ.h"
#include "CpptrajStdio.h"

DataSet_Vector_XYZ::DataSet_Vector_XYZ() :
  DataSet_Vector(VEC_XYZ, TextFormat(TextFormat::DOUBLE, 8, 4, 3))
{}

// DataSet_Vector_XYZ::MemUsageInBytes()
size_t DataSet_Vector_XYZ::MemUsageInBytes() const {
  return internalSize();
}

// DataSet_Vector_XYZ::Allocate
int DataSet_Vector_XYZ::Allocate(SizeArray const& Nin) {
  if (!Nin.empty())
    internalAlloc(Nin);
  return 0;
}

// DataSet_Vector_XYZ::WriteBuffer
void DataSet_Vector_XYZ::WriteBuffer(CpptrajFile& cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Size()) {
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0); // VXYZ
  } else {
    Vec3 const& Vxyz = vectors()[pIn[0]];
    cbuffer.Printf(format_.fmt(), Vxyz[0], Vxyz[1], Vxyz[2]);
  }
}

// DataSet_Vector_XYZ::Append()
int DataSet_Vector_XYZ::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != VECTOR_1D) return 1;
  if (dsIn->Type() == VEC_OXYZ)
    mprintf("Warning: Appending vector set with origins to set without origins; will lose origin info.\n");
  internalAppend( (DataSet_Vector*)dsIn );

  return 0;
}

#ifdef MPI
int DataSet_Vector_XYZ::Sync(size_t total, std::vector<int> const& rank_frames,
                             Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  // TODO: Consolidate to 1 send/recv via arrays?
  if (commIn.Master()) {
    // Resize to accept data from other ranks.
    internalVecArray().resize( total );
    int vidx = rank_frames[0]; // Index on master
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int ridx = 0; ridx != rank_frames[rank]; ridx++, vidx++)
        commIn.SendMaster( internalVecArray()[vidx].Dptr(), 3, rank, MPI_DOUBLE );
    } 
  } else { // Send data to master
    for (unsigned int ridx = 0; ridx != internalVecArray().size(); ++ridx)
      commIn.SendMaster( internalVecArray()[ridx].Dptr(), 3, commIn.Rank(), MPI_DOUBLE );
  }
  return 0;
}
#endif

// -----------------------------------------------------------------------------
void DataSet_Vector_XYZ::reset() {
  internalReset();
}
