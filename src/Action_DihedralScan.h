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
      std::vector<int> checkAtoms;
      double interval;
      int currentVal;
      int maxVal;
      int atom1;
      int atom2;
      int resnum;
      bool isRandom;
    };
    std::vector<DihedralScanType> BB_dihedrals;

    struct ResidueCheckType {
      int checkatom;
      int start;
      int stop;
      int resnum;
    };
    std::vector<ResidueCheckType> ResCheck;

    AtomMask Mask1;
    bool random_angle;
    bool check_for_clashes;
    char *outfilename;
    char *outfmt;
    double interval;
    int max_rotations;
    int max_factor;
    double cutoff;
    double rescutoff;
    int backtrack;
    int increment;     ///< Value in degrees to increment random dihedral by if clash happens
    int max_increment; ///< 360 / increment
    DataSet *number_of_problems;
    CheckStructure checkStructure;

    //int CheckResidues( Frame *, int );
    int CheckResidue( Frame *, DihedralScanType&,int,double*);
  public:
    DihedralScan();

    int init();
    int setup();
    int action();
};
#endif  
