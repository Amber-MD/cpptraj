#ifndef INC_EXEC_PERMUTEDIHEDRALS_H
#define INC_EXEC_PERMUTEDIHEDRALS_H
#include "Exec.h"
#include "Action_CheckStructure.h"
#include "Random.h"
#include "Trajout_Single.h"
#include "DihedralSearch.h"
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

    int GetDihedralIdxs(int*, Topology const&,int, NameType const&,
                         NameType const&, NameType const&);
    int CheckResidue( Frame const&, Topology const&, PermuteDihedralsType const&,int,double*);
    void RandomizeAngles(Frame&, Topology const&);
    void IntervalAngles(Frame const&, Topology const&);

    /// Scan types
    enum ModeType  { RANDOM = 0, INTERVAL };
    /// What kind of scanning will be performed
    ModeType mode_;
    /// Used to search for and hold info for specified dihedrals
    DihedralSearch dihSearch_;
    std::vector<PermuteDihedralsType> BB_dihedrals_;
    /// Hold info for clash check
    struct ResidueCheckType {
      int checkatom;
      int start;
      int stop;
      int resnum;
    };
    std::vector<ResidueCheckType> ResCheck_;

    Range resRange_;
    Trajout_Single outtraj_;
    std::string outfilename_;
    int outframe_;
    double interval_;   ///< interval, value to shift by
    int maxVal_;        ///< Maximum number of times to rotate dihedral
    // 'random' options
    bool check_for_clashes_;
    bool checkAllResidues_;
    int max_rotations_; ///< Max # of random rotations to try, == # of dihedrals
    int max_factor_;    ///< # of times to randomly rotate each dihedral
    double cutoff_;     ///< When checking for clashes, atom cutoff
    double rescutoff_;  ///< When checking for clashes, residue cutoff
    int backtrack_;     ///< When a clash cannot be resolved, # of dihedrals to backtrack
    int increment_;     ///< Value in degrees to increment random dihedral by if clash happens
    int max_increment_; ///< 360 / increment
    // General
    int debug_;
    DataSet* number_of_problems_;
    Action_CheckStructure checkStructure_;
    Random_Number RN_;
};
#endif
