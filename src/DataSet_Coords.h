#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "Topology.h"
#include "DataSet_1D.h"
// NOTE: Should this be a 4D DataSet?
/// Interface to COORDS data sets
class DataSet_Coords : public DataSet_1D {
  public:
    DataSet_Coords() : DataSet_1D(COORDS, 8, 3), numBoxCrd_(0), numVel_(0) {}
    virtual ~DataSet_Coords() {}
    /// Allocate a Frame that can be used to store COORDS 
    Frame AllocateFrame() const;
    /// Add given Frame to this COORDS
    virtual void AddFrame(Frame const&) = 0;
    /// Set COORDS at specified position with Frame
    virtual void SetCRD(int, Frame const&) = 0;
    /// Set given Frame with COORDS at specified position
    virtual void GetFrame(int, Frame&) const = 0;
    /// Set given Frame with COORDS at specified position according to mask
    virtual void GetFrame(int, Frame&, AtomMask const&) const = 0;
    /// Set main topology that will be associated with frames to/from this COORDS
    void SetTopology(Topology const&);
    inline Topology const& Top() const { return top_; }
  protected:
    Topology top_;                     ///< Topology corresponding to coordinates.
    size_t numBoxCrd_;                 ///< Number of box coords (0 or 6).
    size_t numVel_;                    ///< Number of velocities.
};
#endif
