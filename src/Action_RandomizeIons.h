#ifndef INC_ACTION_RANDOMIZEIONS_H
#define INC_ACTION_RANDOMIZEIONS_H
#include "Action.h"
#include "ImagedAction.h"
class Action_RandomizeIons : public Action, ImagedAction {
  public:
    Action_RandomizeIons();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_RandomizeIons(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ImagedAction image_; ///< Imaging routines.
    AtomMask ions_;      ///< Mask of ions to be moved.
    double overlap_;     ///< The minimum allowed distance between ions
    double min_;         ///< The minimum distance to the 'around' mask
    AtomMask around_;    ///< The 'around' mask (region of space for ions to avoid)
    int seed_;           ///< random seed
    int debug_;
    int n_solvent_;
    // TODO: Combine the below 3 into a struct?
    std::vector<int> solventStart_; ///< Solvent molecule start atoms.
    std::vector<int> solventEnd_;   ///< Solvent molecule end atoms.
    std::vector<bool> solvent_;     ///< True if solvent mol begin considered for swap.
};
#endif
