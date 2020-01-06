#include <algorithm> // copy
#include "DataSet_Vector_OXYZ.h"

// CONSTRUCTOR
DataSet_Vector_OXYZ::DataSet_Vector_OXYZ() :
  DataSet_Vector(VEC_OXYZ, TextFormat(TextFormat::DOUBLE, 8, 4, 6))
{}

// DataSet_Vector_OXYZ::MemUsageInBytes()
size_t DataSet_Vector_OXYZ::MemUsageInBytes() const {
  return internalSize() + (origins_.size() * Vec3::DataSize());
}
         
// DataSet_Vector_OXYZ::Allocate()
int DataSet_Vector_OXYZ::Allocate(SizeArray const& Nin) {
  if (!Nin.empty()) {
    internalAlloc(Nin); 
    origins_.reserve( Nin[0] );
  }
  return 0;
}

int DataSet_Vector_OXYZ::MemAlloc(SizeArray const& Nin) {
  if (!Nin.empty()) {
    internalMemalloc(Nin);
    origins_.resize( Nin[0] );
  }
  return 0;
}

void DataSet_Vector_OXYZ::CopyBlock(size_t startIdx, DataSet const* dptrIn, size_t pos, size_t nelts)
{
  DataSet_Vector_OXYZ const& setIn = static_cast<DataSet_Vector_OXYZ const&>( *dptrIn );
  internalCopyBlock(startIdx, setIn, pos, nelts);
  Varray::const_iterator ptr = setIn.origins_.begin() + pos;
  std::copy( ptr, ptr + nelts, origins_.begin() + startIdx );
}

// DataSet_Vector_OXYZ::WriteBuffer()
void DataSet_Vector_OXYZ::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= Size()) {
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); // VXYZ OXYZ
  } else {
    Vec3 const& Vxyz = vectors()[pIn[0]];
    Vec3 const& Oxyz = origins_[pIn[0]];
    cbuffer.Printf(format_.fmt(), Vxyz[0], Vxyz[1], Vxyz[2],
                                  Oxyz[0], Oxyz[1], Oxyz[2]);
  }
}

//  DataSet_Vector_OXYZ::Append()
int DataSet_Vector_OXYZ::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Group() != VECTOR_1D) return 1;
  internalAppend( (DataSet_Vector*)dsIn );

  if (dsIn->Type() == VEC_OXYZ) {
    // incoming set has origins
    Varray const& oIn = ((DataSet_Vector_OXYZ*)dsIn)->origins_;
    size_t oldsize = origins_.size();
    origins_.resize( oldsize + oIn.size() );
    std::copy( oIn.begin(), oIn.end(), origins_.begin() + oldsize );
  } else {
    // no origins in incoming set; use zero
    origins_.resize( Size(), Vec3(0.0) );
  }
  return 0;
}

// DataSet_Vector_OXYZ::reset()
void DataSet_Vector_OXYZ::reset() {
  internalReset(); 
  origins_.clear();
}

#ifdef MPI
int DataSet_Vector_OXYZ::Sync(size_t total, std::vector<int> const& rank_frames,
                              Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  double buf[6];
  // TODO: Consolidate to 1 send/recv via arrays?
  if (commIn.Master()) {
    // Resize to accept data from other ranks.
    internalVecArray().resize( total );
    origins_.resize( total );
    int vidx = rank_frames[0]; // Index on master
    for (int rank = 1; rank < commIn.Size(); rank++) {
      for (int ridx = 0; ridx != rank_frames[rank]; ridx++, vidx++) {
        commIn.SendMaster( buf, 6, rank, MPI_DOUBLE );
        std::copy( buf,   buf+3, internalVecArray()[vidx].Dptr() );
        std::copy( buf+3, buf+6, origins_[vidx].Dptr() );
      }
    }
  } else { // Send data to master
    for (unsigned int ridx = 0; ridx != internalVecArray().size(); ++ridx) {
      std::copy( internalVecArray()[ridx].Dptr(), internalVecArray()[ridx].Dptr()+3, buf   );
      std::copy( origins_[ridx].Dptr(), origins_[ridx].Dptr()+3, buf+3 );
      commIn.SendMaster( buf, 6, commIn.Rank(), MPI_DOUBLE );
    }
  }
  return 0;
}
#endif
