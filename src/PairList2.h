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
    /// PairList Atom
    class AtmType;
    typedef std::vector<AtmType> Aarray;
    /// PairList Cell
    class CellType;
    typedef std::vector<CellType> Carray;

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
    /// \return Number of grid cells.
    int NGridMax()                const { return (int)cells_.size(); }
    /// \return Specified grid cell.
    CellType const& Cell(int idx) const { return cells_[idx];        }
    /// \return Translation vector for given translation index (from TransList()).
    Vec3 const& TransVec(int t)    const { return translateVec_[t];  }
  private:
    /// Determine neighbors and translation vectors for each cell.
    void CalcGridPointers(int,int);
    /// Check grid dimensions using given recip lengths and (re)allocate mem if necessary.
    int SetupGrids(Vec3 const&);
    /// Update the translation vectors based on given unit cell matrix.
    void FillTranslateVec(Matrix_3x3 const&);
    /// Add atoms to grid cells.
    void GridUnitCell(Frame const&, Matrix_3x3 const&, Matrix_3x3 const&, AtomMask const&);
    /// Calculate cell index based on given coords and add to cell
    inline void GridAtom(int, Vec3 const&, Vec3 const&);

    Carray cells_;            ///< Hold all cells in grid
    Vec3 translateVec_[18];   ///< Translate vector array
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
// -----------------------------------------------------------------------------
/** PairList Atom. Holds atom index and wrapped coordinates. */
class PairList2::AtmType {
  public:
    AtmType() : idx_(-1) {}
    /// CONSTRUCTOR - Atom index, wrapped fractional coords, wrapped Cart. coords
    AtmType(int i, Vec3 const& f, Vec3 const& c) :
      imageCoords_(c), fracCoords_(f), idx_(i) {}
    /// \return Wrapped Cart. coords
    Vec3 const& ImageCoords() const { return imageCoords_; }
    /// \return Atom index.
    int Idx()                 const { return idx_;         }
  private:
    Vec3 imageCoords_; ///< Imaged Cartesian coordinates
    Vec3 fracCoords_;  ///< Fractional coordinates
    int idx_;          ///< Atom index
};
// -----------------------------------------------------------------------------
/** PairList Cell. Holds list of atoms and neighbor cells. */
class PairList2::CellType {
  public:
    CellType() {}
    /// Add given atom to this cell
    void AddAtom(AtmType const& a) { atoms_.push_back( a ); }
    /// Clear all atoms in this cell
    void ClearAtoms()              { atoms_.clear();        }
    /// \return Number of atoms in the cell
    unsigned int NatomsInGrid() const { return atoms_.size();  }
    /// \return List of cell neighbors - starts with this cell
    Iarray const& CellList()    const { return neighborPtr_;   }
    /// \return List of translations for neighbors - starts with this cell
    Iarray const& TransList()   const { return neighborTrans_; }
    /// \return Memory used by this cell in bytes
    size_t MemSize() const {
      return ((neighborPtr_.size() + neighborTrans_.size()) * sizeof(int)) +
             (2 * sizeof(Iarray)) + (atoms_.size() * sizeof(AtmType)) + sizeof(Aarray);
    }
    /// Atom list iterator
    typedef Aarray::const_iterator const_iterator;
    /// \return Iterator to beginning of atom list.
    const_iterator begin() const { return atoms_.begin(); }
    /// \return Iterator to end of atom list.
    const_iterator end()   const { return atoms_.end();   }
    // PairList2 is a friend so it can access neighborPtr_ and neighborTrans_
    friend class PairList2;
  private:
    Iarray neighborPtr_;   ///< Indices of neighbor cells "forward" of this cell.
    Iarray neighborTrans_; ///< Index pointing to translate vector for each neighbor cell.
    Aarray atoms_;         ///< Atoms in this cell.
};
#endif
