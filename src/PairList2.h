#ifndef INC_PAIRLIST2_H
#define INC_PAIRLIST2_H
#include "Topology.h"
#include "Timer.h"
/// Class for creating a lists of potential pairing atoms via spatial grid.
/** NOTE: The code in this class is based on that from the SANDER
  *       program of Amber/AmberTools, particularly nonbond_list.F90.
  *       However, it is more memory-hungry since it stores full cell
  *       neighbor lists instead of cell centers.
  */
class PairList2 {
  public:
    typedef std::vector<int> Iarray;
    class Atm;
    typedef std::vector<Atm> Aarray;
    class Cell;
    typedef std::vector<Cell> Carray;

    PairList2();
    /// Initialize pair list with given cutoff, "skin", and debug level.
    int InitPairList(double,double,int);
    /// Setup pair list grid cells based on given box and vector of recip lengths.
    int SetupPairList(Box::BoxType, Vec3 const&);
    /// Create pair list from Frame, unit cell and recip matrices, and mask.
    int CreatePairList(Frame const&, Matrix_3x3 const&, Matrix_3x3 const&, AtomMask const&);
    /// Print timing info.
    void Timing(double) const;
    /// Print memory usage.
    void PrintMemory() const;

  private:
    /// Determine neighbors and translation vectors for each cell.
    void CalcGridPointers(int,int);
    /// Check grid dimensions using given recip lengths and (re)allocate mem if necessary.
    int SetupGrids(Vec3 const&);


    Carray cells_;            ///< Hold all cells in grid

    double cutList_;          ///< Direct space cutoff plus non-bond "skin"
    int debug_;
    int nGridX_;              ///< Number of grid cells in X direction.
    int nGridY_;              ///< Number of grid cells in Y direction.
    int nGridZ_;              ///< Number of grid cells in Z direction.
    int nGridX_0_;            ///< Previous number of cells in X direction
    int nGridY_0_;            ///< Previous number of cells in Y direction
    int nGridZ_0_;            ///< Previous number of cells in Z direction

    static const int cellOffset_; ///< Number of cells in forward direction to check

    Timer t_map_;
    Timer t_gridpointers_;
    Timer t_total_;
};
/// PairList Atom
class PairList2::Atm {
  public:
    Atm() {}
  private:
    Vec3 imageCoords_; ///< Imaged Cartesian coordinates
    Vec3 fracCoords_;  ///< Fractional coordinates
    int idx_;          ///< Atom index
};
/// PairList Cell
class PairList2::Cell {
  public:
    Cell() {}
  private:
    Iarray neighborPtr_;   ///< Indices of neighbor cells "forward" of this cell.
    Iarray neighborTrans_; ///< Index pointing to translate vector for each neighbor cell.
    Aarray atoms_;         ///< Atoms in this cell.
};
#endif
