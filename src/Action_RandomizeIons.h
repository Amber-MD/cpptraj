#ifndef INC_ACTION_RANDOMIZEIONS_H
#define INC_ACTION_RANDOMIZEIONS_H
#include "Action.h"
#include "ImagedAction.h"
#include "Random.h"
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

    ImagedAction image_; ///< Imaging routines.
    Random_Number RN_;   ///< Random number generator.
    AtomMask ions_;      ///< Mask of ions to be moved.
    AtomMask around_;    ///< The 'around' mask (region of space for ions to avoid)
    double overlap_;     ///< The minimum allowed distance between ions
    double min_;         ///< The minimum distance to the 'around' mask
    int n_solvent_;      ///< Total number of solvent molecules.
    int debug_;
    // TODO: Combine the below 3 into a struct?
    std::vector<int> solventStart_; ///< Solvent molecule start atoms.
    std::vector<int> solventEnd_;   ///< Solvent molecule end atoms.
    std::vector<bool> solvent_;     ///< True if solvent mol being considered for swap.
};
#endif
