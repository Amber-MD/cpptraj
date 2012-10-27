#ifndef INC_ACTION_DIHEDRALSCAN_H
#define INC_ACTION_DIHEDRALSCAN_H
#include "Action.h"
#include "Action_CheckStructure.h"
#include "Random.h"
// Class: Action_DihedralScan
/// Action to rotate dihedrals randomly or in intervals 
class Action_DihedralScan: public Action {
  public:
    Action_DihedralScan();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_DihedralScan(); }
    static void Help();

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
  private:
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
    std::string outfilename_;
    std::string outfmt_;
    double interval;
    int max_rotations;
    int max_factor;
    double cutoff;
    double rescutoff;
    int backtrack;
    int increment;     ///< Value in degrees to increment random dihedral by if clash happens
    int max_increment; ///< 360 / increment
    int debug_;
    Topology* CurrentParm_;
    DataSet *number_of_problems;
    Action_CheckStructure checkStructure;
    Random_Number RN_;

    //int CheckResidues( Frame *, int );
    int CheckResidue( Frame *, DihedralScanType&,int,double*);
};
#endif  
