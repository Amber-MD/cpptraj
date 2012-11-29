#ifndef INC_ACTION_DIHEDRALSCAN_H
#define INC_ACTION_DIHEDRALSCAN_H
#include "Action.h"
#include "Action_CheckStructure.h"
#include "Random.h"
#include "Trajout.h"
// Class: Action_DihedralScan
/// Action to rotate dihedrals randomly or in intervals 
class Action_DihedralScan: public Action {
  public:
    Action_DihedralScan();
    ~Action_DihedralScan();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_DihedralScan(); }
    static void Help();

    void Print() {}
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    /// What kind of scanning will be performed
    enum ModeType { RANDOM = 0, INTERVAL, IMPOSE };
    ModeType mode_;
    /// Hold info for a dihedral
    struct DihedralScanType {
      AtomMask Rmask;
      std::vector<int> checkAtoms;
      double interval;
      int currentVal;
      int maxVal;
      int atom1;
      int atom2;
      int resnum;
    };
    std::vector<DihedralScanType> BB_dihedrals_;

    struct ResidueCheckType {
      int checkatom;
      int start;
      int stop;
      int resnum;
    };
    std::vector<ResidueCheckType> ResCheck_;

    AtomMask Mask1_;
    bool check_for_clashes_;
    Trajout outtraj_;
    std::string outfilename_;
    int outframe_;
    double interval_;
    int max_rotations_;
    int max_factor_;
    double cutoff_;
    double rescutoff_;
    int backtrack_;
    int increment_;     ///< Value in degrees to increment random dihedral by if clash happens
    int max_increment_; ///< 360 / increment
    int debug_;
    Topology* CurrentParm_;
    DataSet *number_of_problems_;
    Action_CheckStructure checkStructure_;
    Random_Number RN_;

    //int CheckResidues( Frame *, int );
    int CheckResidue( Frame const&, DihedralScanType&,int,double*);
    void IntervalAngles(Frame&);
    void RandomizeAngles(Frame&);
};
#endif  
