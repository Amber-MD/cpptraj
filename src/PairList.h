#ifndef INC_PAIRLIST_H
#define INC_PAIRLIST_H
#include "Topology.h"
#include "Timer.h"
/// Class for creating a lists of potential pairing atoms via spatial grid cutoff.
/** NOTE: The code in this class is largely based on that from the SANDER
  *       program of Amber/AmberTools, particularly nonbond_list.F90.
  */
class PairList {
  public:
    PairList();
    int InitPairList();
    void CreatePairList(Frame const&, AtomMask const&);
  private:
    int Fill_CellNeighbor();
    void MapCoords(Frame const&, Matrix_3x3 const&,Matrix_3x3 const&, AtomMask const&);
    void FillTranslateVec(Matrix_3x3 const&);

    typedef std::vector<Vec3> Varray;

    static const int cellOffset_; ///< Number of cells in forward direction to check

    Vec3 translateVec_[18];   ///< Translate vector array
    int cellNeighbor_[7][10]; ///< Cell neighbor index array based on cellOffset_.
    Varray Frac_;             ///< Hold fractional coords back in primary cell.
    Varray Image_;            ///< Hold Cartesian coords back in primary cell.
    Timer t_map_;
};
#endif
