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
    int InitPairList(double,double);
    int CreatePairList(Frame const&, AtomMask const&);
  private:
    int Fill_CellNeighbor();
    void MapCoords(Frame const&, Matrix_3x3 const&,Matrix_3x3 const&, AtomMask const&);
    void FillTranslateVec(Matrix_3x3 const&);
    int SetupGrids(Vec3 const&);
    void GridUnitCell();

    typedef std::vector<Vec3> Varray;
    typedef std::vector<int> Iarray;

    static const int cellOffset_; ///< Number of cells in forward direction to check

    Vec3 translateVec_[18];   ///< Translate vector array
    int cellNeighbor_[7][10]; ///< Cell neighbor index array based on cellOffset_.
    Varray Frac_;             ///< Hold fractional coords back in primary cell.
    Varray Image_;            ///< Hold Cartesian coords back in primary cell.
    Iarray nLoGrid_;
    Iarray nHiGrid_;
    Iarray MyGrids_;
    Iarray nAtomsInGrid_;     ///< Number of atoms in each grid cell.
    Iarray idxOffset_;        ///< Offset of starting atom in each grid cell.
    Iarray atomCell_;         ///< Grid cell index for each atom.
    double cutList_;          ///< Direct space cutoff plus non-bond "skin"
    int nGridX_;              ///< Number of grid cells in X direction.
    int nGridY_;              ///< Number of grid cells in Y direction.
    int nGridZ_;              ///< Number of grid cells in Z direction.
    int nGridMax_;            ///< Total number of grid cells.
    Timer t_map_;
};
#endif
