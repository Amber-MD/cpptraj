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
    DataSet *surf;
    AtomMask Mask1;
    AtomMask atomi_neighborMask;
    AtomMask atomi_noNeighborMask;
    AtomMask atomj_neighborMask;

  public:
    Surf();

    int init();
    int setup();
    int action();
};
#endif
