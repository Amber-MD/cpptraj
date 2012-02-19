#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include <vector>
#include "Action.h"
// Class: CheckStructure 
/// Action to check bond lengths and bad overlaps between non-bonded atoms 
class CheckStructure: public Action {
    /// Used to cache bond parameters
    struct bond_list {
      double req;
      int atom1;
      int atom2;
      int param;
    };
    std::vector<bond_list> bondL;
    // Sort first by atom1, then by atom2
    struct bond_list_cmp {
      inline bool operator()(bond_list first, bond_list second) const {
        if (first.atom1 < second.atom1) {
          return true;
        } else if (first.atom1 == second.atom1) {
          if (first.atom2 < second.atom2) return true;
        } 
        return false;
      }
    };

    AtomMask Mask1;
    double bondoffset;
    double nonbondcut2;
    CpptrajFile outfile;
  public:
    CheckStructure();
    ~CheckStructure();

    int init();
    int setup();
    int action();

    void SeparateInit(double, double, int);
    int SeparateAction(Frame *);
};
#endif  
