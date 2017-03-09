#ifndef INC_PAIRLIST_H
#define INC_PAIRLIST_H
#include "Topology.h"
#include "Timer.h"
/// Class for creating a lists of potential pairing atoms via spatial grid.
/** NOTE: The code in this class is based on that from the SANDER
  *       program of Amber/AmberTools, particularly nonbond_list.F90.
  *       However, it is more memory-hungry since it stores full cell
  *       neighbor lists instead of cell centers.
  */
class PairList {
  public:
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Iarray2D;
    typedef std::vector<Vec3> Varray;
    PairList();
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
    /// \return Total number of grid cells.
    int NGridMax()                 const { return nGridMax_;         }
    /// \return Array of neighbor grid cell indices for specified grid cell.
    Iarray const& Cell(int i)      const { return neighborPtr_[i];   }
    /// \return Array of indices into TransVec for each neighbor grid cell in Cell(i).
    Iarray const& Trans(int i)     const { return neighborTrans_[i]; }
    /// \return Starting index into AtomGridIdx for given grid cell.
    int IdxOffset(int i)           const { return idxOffset_[i];     }
    /// \return Number of atoms in given grid cell.
    int NatomsInGrid(int i)        const { return nAtomsInGrid_[i];  }
    /// \return Atom index for given index (related to IdxOffset)
    int AtomGridIdx(int a)         const { return atomGridIdx_[a];   }
    /// \return Imaged coordinates for given atom index.
    Vec3 const& ImageCoords(int n) const { return Image_[n];         }
    /// \return Translation vector for given translation index (from Trans()).
    Vec3 const& TransVec(int t)    const { return translateVec_[t];  }

    Varray const& ImageCoords()    const { return Image_; }
    Varray const& FracCoords()     const { return Frac_;  }
  private:
//    int Fill_CellNeighbor();
    /// Determine neighbors and translation vectors for each cell.
    void CalcGridPointers(int,int);
    /// Check grid dimensions using given recip lengths and (re)allocate mem if necessary.
    int SetupGrids(Vec3 const&);
    /// Convert selected coords to wrapped fraction coords and save wrapped Cartesian coords.
    void MapCoords(Frame const&, Matrix_3x3 const&,Matrix_3x3 const&, AtomMask const&);
    /// Update the translation vectors based on given unit cell matrix.
    void FillTranslateVec(Matrix_3x3 const&);
    /// Assign fraction coords to grid cells.
    void GridUnitCell();

//    typedef std::vector<bool> Barray;

    static const int cellOffset_; ///< Number of cells in forward direction to check

    Vec3 translateVec_[18];   ///< Translate vector array
//    int cellNeighbor_[7][10]; ///< Cell neighbor index array based on cellOffset_.
    Iarray2D neighborPtr_;    ///< For each grid cell, hold indexes into neighbor cells.
    Iarray2D neighborTrans_;
    Varray Frac_;             ///< Hold fractional coords back in primary cell.
    Varray Image_;            ///< Hold Cartesian coords back in primary cell.
//    Iarray nLoGrid_;
//    Iarray nHiGrid_;
    Iarray nAtomsInGrid_;     ///< Number of atoms in each grid cell.
    Iarray idxOffset_;        ///< Offset of starting atom in each grid cell.
    Iarray atomCell_;         ///< Grid cell index for each atom.
    Iarray atomGridIdx_;      ///< List of atoms sorted by grid cell.
//    Barray myGrids_;          ///< True if I am responsible for this grid cell.
    double cutList_;          ///< Direct space cutoff plus non-bond "skin"
    int debug_;
    int NP_;                  ///< Total number of elements in neighborPtr_/neighborTrans_
    int nGridX_;              ///< Number of grid cells in X direction.
    int nGridY_;              ///< Number of grid cells in Y direction.
    int nGridZ_;              ///< Number of grid cells in Z direction.
    int nGridMax_;            ///< Total number of grid cells.
    int nGridX_0_;            ///< Previous number of cells in X direction
    int nGridY_0_;            ///< Previous number of cells in Y direction
    int nGridZ_0_;            ///< Previous number of cells in Z direction
//    int maxNptrs_;            ///< Max number of neighbor pointers.
    Timer t_map_;
    Timer t_gridpointers_;
    Timer t_total_;
};
#endif
