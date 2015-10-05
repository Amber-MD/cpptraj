#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "DataSet.h"
#include "Topology.h"
/// Interface to COORDS data sets
class DataSet_Coords : public DataSet {
  public:
    DataSet_Coords() {}
    DataSet_Coords(DataType t) : DataSet(t, COORDINATES, TextFormat(), 1) {}
    virtual ~DataSet_Coords() {}
    // -------------------------------------------
    // NOTE: Disabled for all COORDS style DataSets
    void WriteBuffer(CpptrajFile&, SizeArray const&) const {}
    int Append(DataSet*) { return 1; }
    // -------------------------------------------
    /// Add given Frame to this COORDS
    virtual void AddFrame(Frame const&) = 0;
    /// Set COORDS at specified position with Frame
    virtual void SetCRD(int, Frame const&) = 0;
    /// Set given Frame with COORDS at specified position
    virtual void GetFrame(int, Frame&) = 0;
    /// Set given Frame with COORDS at specified position according to mask
    virtual void GetFrame(int, Frame&, AtomMask const&) = 0;
    /// Set topology and coordinate information associated with this COORDS set.
    virtual int CoordsSetup(Topology const&, CoordinateInfo const&) = 0;
    // -------------------------------------------
    /// Allocate a Frame that can be used to store COORDS 
    Frame AllocateFrame() const;
    /// \return Topology associated with these COORDS.
    inline Topology const& Top() const { return top_; }
    /// \return Pointer to Topology associated with these COORDS.
    inline Topology* TopPtr() { return &top_; }
    /// \return CoordinateInfo associated with these COORDS
    inline CoordinateInfo const& CoordsInfo() const { return cInfo_; }
  protected:
    void CommonInfo() const;
    Topology top_;         ///< Topology corresponding to coordinates.
    CoordinateInfo cInfo_; ///< Describes coordinate Frame
};
#endif
