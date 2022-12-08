#ifndef INC_ACTION_BOUNDS_H
#define INC_ACTION_BOUNDS_H
#include "Action.h"
#include <cstddef> // size_t
/// Report the min/max XYZ values for atoms in mask.
class Action_Bounds : public Action {
  public:
    Action_Bounds();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Bounds(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();

    int calc_center_and_bins();

    AtomMask mask_;        ///< Mask of atoms to obtain the bounds of.
    CpptrajFile* outfile_; ///< File to print bounds to.
    double max_[3];        ///< Minimum extent of any atom in mask.
    double min_[3];        ///< Maximum extent of any atom in mask.
    Vec3 dxyz_;            ///< Grid spacing (if creating a grid from bounds).
    size_t nxyz_[3];       ///< # grid bins in each dimension based on spacing.
    Vec3 center_;          ///< Grid center based on min/max.
    int offset_;           ///< # bins offset for grid.
    bool gridIsSetup_;     ///< True if grid has already been set up for bounds.
    DataSet* grid_;        ///< Hold grid created from bounds.
    DataSet* ds_xmin_;
    DataSet* ds_ymin_;
    DataSet* ds_zmin_;
    DataSet* ds_xmax_;
    DataSet* ds_ymax_;
    DataSet* ds_zmax_;
};
#endif
