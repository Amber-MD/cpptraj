#ifndef INC_PAIRLIST_H
#define INC_PAIRLIST_H
#include "Topology.h"
/// Class for creating a lists of potential pairing atoms via spatial grid cutoff.
/** NOTE: The code in this class is largely based on that from the SANDER
  *       program of Amber/AmberTools, particularly nonbond_list.F90.
  */
class PairList {
  public:
    PairList();
  private:
    int Fill_CellNeighbor();

    typedef std::vector<Vec3> Varray;

    static const int cellOffset_; ///< Number of cells in forward direction to check

    Vec3 translateVec_[18];   ///< Translate vector array
    int cellNeighbor_[7][10]; ///< Cell neighbor index array based on cellOffset_.
};
#endif
