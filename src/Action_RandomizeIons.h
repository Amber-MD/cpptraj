#ifndef INC_ACTION_RANDOMIZEIONS_H
#define INC_ACTION_RANDOMIZEIONS_H
#include "Action.h"
class Action_RandomizeIons : public Action {
  public:
    Action_RandomizeIons();
  private:
    int init();
    int setup();
    int action();

    AtomMask ions_;    ///< the list of ions to be moved.
    double overlap_;   ///< darg1: the minimum distance between ions
    double min_;       ///< darg2: the minimum distance to the around mask
    AtomMask around_;  ///< carg1: the around mask (region of space to avoid)
    std::string aroundmask_; ///< empty if no around mask specified.
    int seed_;         ///< iarg2: random seed
    // TODO: Combine the below 3 into a struct?
    /// Hold solvent molecule start atoms.
    std::vector<int> solventStart_;
    /// Hold solvent molecule end atoms.
    std::vector<int> solventEnd_;
    /// True is solvent mol being considered for swap.
    std::vector<bool> solvent_;
};
#endif
