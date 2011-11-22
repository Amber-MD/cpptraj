#ifndef INC_ACTION_DIHEDRALSCAN_H
#define INC_ACTION_DIHEDRALSCAN_H
#include <vector>
#include "Action.h"
#include "Action_CheckStructure.h"
// Class: DihedralScan
/// Action to rotate dihedrals randomly or in intervals 
class DihedralScan: public Action {
    struct DihedralScanType {
      AtomMask Rmask;
      double interval;
      int currentVal;
      int maxVal;
      int atom1;
      int atom2;
      bool isRandom;
    };
    std::vector<DihedralScanType> BB_dihedrals;

    AtomMask Mask1;
    bool random_angle;
    char *outfilename;
    char *outfmt;
    double interval;
    CheckStructure checkStructure;
  public:
    DihedralScan();
    ~DihedralScan();

    int init();
    int setup();
    int action();
};
#endif  
