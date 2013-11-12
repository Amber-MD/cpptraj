#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "DataSet_1D.h"
#include "Topology.h"
// NOTE: Should this be a 4D DataSet?
class DataSet_Coords : public DataSet_1D {
  public:
    DataSet_Coords();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords();    }
    // ----- DataSet functions -------------------
    size_t Size() const { return coords_.size(); }
    int Sync()          { return 1;              }
    void Info() const;
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* )              { return;     }
    double Dval(size_t)                    const { return 0.0; }
    double Xcrd(size_t idx)  const { return Dim(0).Coord(idx); }
    void WriteBuffer(CpptrajFile&, size_t) const { return;     }
    // -------------------------------------------
    /// Allocate a frame that can be used to store COORDS
    Frame AllocateFrame() const {
      Frame f;
      f.SetupFrameV( top_.Atoms(), (numVel_ > 0), 0 );
      return f;
    }
    /// Add a frame.
    inline void AddFrame(Frame const& fIn) { 
      coords_.push_back( fIn.ConvertToCRD(numVel_, numBoxCrd_) ); 
    }
    /// Get a frame at position.
    inline void GetFrame(int idx, Frame& fIn) const { 
      fIn.SetFromCRD( coords_[idx] ); 
    }
    /// Get a frame at position corresponding to mask.
    inline void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) const {
      fIn.SetFromCRD( coords_[idx], mIn );
    }
    /// Set CRD at position with frame.
    inline void SetCRD(int idx, Frame const& fIn) {
      coords_[idx] = fIn.ConvertToCRD(numVel_, numBoxCrd_);
    }
    /// Set topology and number of box coords.
    void SetTopology(Topology const&);
    /// \return Topology corresponding to coords.
    inline Topology const& Top() const { return top_; }
  private:
    typedef std::vector<CRDtype> CRDarray;
    CRDarray coords_;                  ///< Array of coordinate frames.
    Topology top_;                     ///< Topology corresponding to coordinates.
    size_t numBoxCrd_;                 ///< Number of box coords (0 or 6).
    size_t numVel_;                    ///< Number of velocities.
};
#endif
