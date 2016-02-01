#ifndef INC_EXEC_PERMUTEDIHEDRALS_H
#define INC_EXEC_PERMUTEDIHEDRALS_H
#include "Exec.h"
#include "Action_CheckStructure.h"
#include "Random.h"
#include "Trajout_Single.h"
#include "DataSet_Coords_CRD.h"
// NOTE: Formerly Action_PermuteDihedrals
/// Rotate dihedrals randomly or in intervals
class Exec_PermuteDihedrals : public Exec {
  public:
    Exec_PermuteDihedrals();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PermuteDihedrals(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    /// Hold additional info for a dihedral
    struct PermuteDihedralsType {
      AtomMask Rmask;              ///< Mask of atoms to hold fixed during rotation
      std::vector<int> checkAtoms; ///< Atoms in same residue that should be checked for clashes
      int atom0;
      int atom1;
      int atom2;
      int atom3;
      int resnum;
    };
    /// Rotate dihedrals by specified interval until original structure is reached
    void IntervalAngles(Frame const&, Topology const&, double);
    /// Used to check for clashes when performing random rotations
    int CheckResidue( Frame const&, Topology const&, PermuteDihedralsType const&,int,double&);
    /// Randomly rotate dihedrals
    void RandomizeAngles(Frame&, Topology const&);

    /// Permute types 
    enum ModeType  { RANDOM = 0, INTERVAL };
    ModeType mode_;            ///< What kind of dihedral permutations will be performed 
    std::vector<PermuteDihedralsType> BB_dihedrals_;
    /// Hold info for clash check
    struct ResidueCheckType {
      int checkatom;
      int start;
      int stop;
      int resnum;
    };
    std::vector<ResidueCheckType> ResCheck_;
    // General
    int debug_;
    Trajout_Single outtraj_;     ///< Output trajectory
    int outframe_;               ///< Output trajectory frame count
    DataSet_Coords_CRD* crdout_; ///< Output COORDS set
    // 'random' options
    bool check_for_clashes_;
    bool checkAllResidues_;
    int max_factor_;    ///< # of times to randomly rotate each dihedral
    double cutoff_;     ///< When checking for clashes, atom cutoff
    double rescutoff_;  ///< When checking for clashes, residue cutoff
    int backtrack_;     ///< When a clash cannot be resolved, # of dihedrals to backtrack
    int increment_;     ///< Value in degrees to increment random dihedral by if clash happens
    int max_increment_; ///< 360 / increment
    DataSet* number_of_problems_;
    Action_CheckStructure checkStructure_;
    Random_Number RN_;
};
#endif
