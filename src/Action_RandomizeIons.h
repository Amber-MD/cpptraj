#ifndef INC_ACTION_RANDOMIZEIONS_H
#define INC_ACTION_RANDOMIZEIONS_H
#include "Action.h"
#include "ImagedAction.h"
#include "Random.h"
/// Used to randomize ion positions by swapping with solvent molecules.
class Action_RandomizeIons : public Action {
  public:
    Action_RandomizeIons();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_RandomizeIons(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ImagedAction image_;         ///< Imaging routines.
    Random_Number RN_;           ///< Random number generator.
    AtomMask ions_;              ///< Mask of ions to be moved.
    AtomMask around_;            ///< The 'around' mask (region of space for ions to avoid)
    double overlap_;             ///< The minimum allowed distance between ions
    double min_;                 ///< The minimum distance to the 'around' mask
    int debug_;
    std::vector<Unit> solvMols_; ///< Hold all solvent molecules
    std::vector<bool> solvent_;  ///< True if solvent mol being considered for swap.
};
#endif
