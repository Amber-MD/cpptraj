#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "DataSet.h"
#include "Topology.h"
class DataSet_Coords : public DataSet {
  public:
    DataSet_Coords();
    // DataSet routines
    int Allocate(int);
    int Size() { return (int)coords_.size(); }
    void Info();
    /// Add a frame.
    void AddFrame(Frame const& fIn) { 
      coords_.push_back( fIn.ConvertToCRD(numBoxCrd_) ); 
    }
    /// Get a frame at position.
    void GetFrame(int idx, Frame& fIn) const { 
      fIn.SetFromCRD( coords_[idx], numBoxCrd_ ); 
    }
    /// Get a frame at position corresponding to mask.
    void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) const {
      fIn.SetFromCRD( coords_[idx], numBoxCrd_, mIn );
    }
    /// Set CRD at position with frame.
    void SetCRD(int idx, Frame const& fIn) {
      coords_[idx] = fIn.ConvertToCRD(numBoxCrd_);
    }
    /// Set topology and number of box coords.
    void SetTopology(Topology const&);
    /// \return Topology corresponding to coords.
    Topology const& Top()                     const { return top_;          }
    /// \return CRD at position.
    const Frame::CRDtype& operator[](int idx) const { return coords_[idx];  }
  private:
    typedef std::vector<Frame::CRDtype> CRDarray;
    CRDarray coords_;                  ///< Array of coordinate frames.
    Topology top_;                     ///< Topology corresponding to coordinates.
    size_t numBoxCrd_;                 ///< Number of box coords (0 or 6).
};
#endif
