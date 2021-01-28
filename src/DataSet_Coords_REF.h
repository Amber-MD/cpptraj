#ifndef INC_DATASET_COORDS_REF_H
#define INC_DATASET_COORDS_REF_H
#include "DataSet_Coords.h"
// Forward declarations
class ArgList;
/// Store a single reference frame in double precision.
class DataSet_Coords_REF : public DataSet_Coords {
  public:
    DataSet_Coords_REF() : DataSet_Coords(REF_FRAME) {}
    static DataSet* Alloc() { return (DataSet*) new DataSet_Coords_REF(); }
    // ----- DataSet functions -------------------
    // NOTE: Technically a 1D data set so return 1 if not empty.
    size_t Size() const { if (!frame_.empty()) return 1; else return 0; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info() const;
    void Add( size_t, const void* )              { return;     }
    // Size is only ever 1, no need to allocate.
    int Allocate(SizeArray const&)               { return 0;   }
    size_t MemUsageInBytes() const { return frame_.DataSize(); }
    // ----- DataSet_Coords functions ------------
    /// Add a frame.
    void AddFrame(Frame const&);
    /// Get a frame at position.
    void GetFrame(int idx, Frame& fIn) { fIn.SetFrame( frame_ ); }
    /// Get a frame at position corresponding to mask.
    void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) {
      fIn.SetFrame(frame_, mIn);
    }
    /// Set CRD at position with frame.
    void SetCRD(int, Frame const&);
    /// Set Topology and coordinate info
    int CoordsSetup(Topology const&, CoordinateInfo const&);
    // -------------------------------------------
    /// Set up reference frame from file.
    int LoadRefFromFile(FileName const&, Topology const&, int);
    /// Set up reference frame from file.
    int LoadRefFromFile(FileName const&, std::string const&, Topology const&, ArgList&, int);
    /// Set up reference frame from COORDS DataSet.
    int SetRefFromCoords(DataSet_Coords*, std::string const&, int);
    int StripRef(std::string const&);
    int StripRef(AtomMask const&);
    Frame const& RefFrame()         const { return frame_; }
  private:
    Frame frame_;       ///< Reference coords.
};
#endif
