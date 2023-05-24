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

    Varray torPositions_;  ///< Current positions of each entity in the toroidal scheme.
    Varray prevPositions_; ///< Previous Cartesian positions of each entity.
    Marray entities_;      ///< Masks selecting each entity to track diffusion of.
    CharMask mask1_;       ///< Mask selecting entities.
    bool useMass_;         ///< Control center of mass vs geometric center
    DataSet* avg_x_;  ///< Hold average diffusion in X direction each frame
    DataSet* avg_y_;  ///< Hold average diffusion in Y direction each frame
    DataSet* avg_z_;  ///< Hold average diffusion in Z direction each frame
    DataSet* avg_r_;  ///< Hold average MSD each frame
    DataSet* avg_a_;  ///< Hold average distance each frame
};
#endif
