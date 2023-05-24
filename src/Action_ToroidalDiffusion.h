#ifndef INC_ACTION_TOROIDALDIFFUSION_H
#define INC_ACTION_TOROIDALDIFFUSION_H
#include "Action.h"
#include "CharMask.h"
/// Implements the Toroidal-view-preserving scheme of Hummer et al. for calculating diffusion
class Action_ToroidalDiffusion : public Action {
  public:
    Action_ToroidalDiffusion();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_ToroidalDiffusion(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef std::vector<Vec3> Varray;
    typedef std::vector<AtomMask> Marray;

    Marray setup_entities(Topology const&) const;

    Varray torPositions_; ///< Current positions of each entity in the toroidal scheme.
    Marray entities_;     ///< Masks selecting each entity to track diffusion of.
    CharMask mask1_;      ///< Mask selecting entities.
    bool useMass_;        ///< Control center of mass vs geometric center
};
#endif
