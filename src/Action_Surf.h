#ifndef INC_ACTION_SURF_H
#define INC_ACTION_SURF_H
#include "Action.h"
#include <vector>
// Class: Surf
/// Calculate LCPO surface area.
/** LCPO method from:
  * -  J. Weiser, P.S. Shenkin, and W.C. Still,
  *    "Approximate atomic surfaces from linear combinations of pairwise
  *    overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  */
class Surf: public Action {
  public:
    Surf();

    int init();
    int setup();
    int action();
  private:
    DataSet *surf;
    AtomMask Mask1;
    AtomMask atomi_neighborMask;
    AtomMask atomi_noNeighborMask;
    AtomMask atomj_neighborMask;
    /// Contain data for an atoms LCPO SA calc
    // TODO: Rework VDW storage
    struct SurfInfo {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    /// Contain LCPO data for all atoms in atomi_neighborMask
    std::vector<SurfInfo> SurfaceInfo_neighbor;
    /// Contain LCPO data for all atoms in atomi_noNeighborMask
    std::vector<SurfInfo> SurfaceInfo_noNeighbor;
    /// Contain vdw radii for all atoms
    std::vector<double> VDW;

    void AssignLCPO(SurfInfo *, double, double, double, double , double );
    void SetAtomLCPO(int,const Atom &, SurfInfo*);
};
#endif
