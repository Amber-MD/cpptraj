#ifndef INC_ACTION_SURF_H
#define INC_ACTION_SURF_H
#include "Action.h"
// Class: Action_Surf
/// Calculate LCPO surface area.
/** LCPO method from:
  * -  J. Weiser, P.S. Shenkin, and W.C. Still,
  *    "Approximate atomic surfaces from linear combinations of pairwise
  *    overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  */
class Action_Surf: public Action {
  public:
    Action_Surf();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Surf(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    /// Struct for storing vdW and LCPO SA params for an atom
    struct SurfInfo {
      double vdwradii;
      double P1;
      double P2;
      double P3;
      double P4;
    };
    /// Assign LCPO vdW radius and parameters 1-4 to given SurfInfo.
    void AssignLCPO(SurfInfo*, double, double, double, double, double);
    /// Assign LCPO vdw and parameters to SurfInfo for specified atom.
    void SetAtomLCPO(Topology const&,int, SurfInfo*);

    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;
    typedef std::vector<SurfInfo> Parray;

    DataSet* surf_;       ///< Hold LCPO surface area data
    AtomMask Mask1_;      ///< Atoms to calculate SA for.
    AtomMask SoluteMask_; ///< Used to select solute atoms.
    Iarray HeavyAtoms_;   ///< Solute atoms with vdW > neighborCut (need vdW radii only)
    Darray VDW_;          ///< Hold vdW radii for HeavyAtoms_
    Iarray SA_Atoms_;     ///< Selected solute atoms with vdW > neighborCut (need vdW + SA params)
    Parray Params_;       ///< Hold vdW + SA params for atoms in SA_Atoms_ 
    double neighborCut_;    ///< Atoms with vdW > this have neighbors.
    double noNeighborTerm_; ///< SA contribution from atoms with no neighbors.
    double offset_;         ///< vdW offset; Amber default is 1.4 Ang.
    /// Hold indices of atoms that are neighbors to the current atom.
#   ifdef _OPENMP
    std::vector<Iarray> Ineighbor_;
#   else
    Iarray Ineighbor_;
#   endif
    /// Hold distances from current atom to neighbors.
#   ifdef _OPENMP
    std::vector<Darray> DIJ_;
#   else
    Darray DIJ_;
#   endif
};
#endif
