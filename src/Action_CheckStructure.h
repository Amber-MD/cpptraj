#ifndef INC_ACTION_CHECKSTRUCTURE_H
#define INC_ACTION_CHECKSTRUCTURE_H
#include <vector>
#include "Action.h"
// Class: CheckStructure 
/// Action to check bond lengths and bad overlaps between non-bonded atoms 
class CheckStructure: public Action {
    std::vector<double> req_array; ///< Hold req for bonded atom pairs, 0 for non bonded
    bool noimage;
    int imageType; 
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
};
#endif  
