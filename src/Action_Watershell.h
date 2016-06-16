#ifndef INC_ACTION_WATERSHELL_H
#define INC_ACTION_WATERSHELL_H
#include "Action.h"
#include "ImagedAction.h"
/// Calculate number of solvent residues in 1st/2nd solvation shell.
class Action_Watershell : public Action {
  public:
    Action_Watershell();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Watershell(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() { return; }

    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    ImagedAction image_;    ///< Hold imaging routines.
    AtomMask soluteMask_;   ///< Selected solute atoms.
    AtomMask solventMask_;  ///< Selected solvent atoms.
    double lowerCutoff_;    ///< Solvent below this is in the first shell.
    double upperCutoff_;    ///< Solvent below this is in the second shell.
    Topology* CurrentParm_; ///< Used to get molecule number for each solvent atom.
    DataSet* lower_;        ///< Number of solvent in first shell.
    DataSet* upper_;        ///< Number of solvent in second shell.
    Darray soluteCoords_;   ///< Hold selected solute coords.
#   ifdef _OPENMP
    /// Shell status for solvent for each OpenMP thread.
    std::vector<Iarray> shellStatus_thread_; ///< Shell status for solvent for each OpenMP thread.
#   else
    Iarray shellStatus_;    ///< Solvent shell status for each solvent (none, first, second)
#   endif
};
#endif
