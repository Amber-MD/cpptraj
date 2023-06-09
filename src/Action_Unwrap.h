#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
#include "ImageTypes.h"
namespace Image {
  class List;
}
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();
    ~Action_Unwrap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unwrap(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum UnwrapType { FRAC = 0, TOR };

    Image::List* imageList_;     ///< List of entities to unwrap
    Image::Mode imageMode_;      ///< Specifiy unwrapping by atom, res, or mol
    std::string maskExpression_; ///< Expression selecting entities to unwrap
    Frame RefFrame_;             ///< Reference frame, updated each DoAction
    Topology* RefParm_;          ///< Reference topology
    bool center_;                ///< If true, determine distances to centers
    bool refNeedsCalc_;          ///< If true, need to calc frac. coords of ref
    Unit allEntities_;           ///< Hold atoms to copy from target to reference
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
    typedef std::vector<Vec3> Varray;
    Varray fracCoords_;   ///< Hold fractional coords for previous frame
    DataSet* avgucell_;   ///< Hold average unit cell parameters for removing box fluctuations
    Box avgbox_;          ///< Hold average box for removing box fluctuations
    Varray torPositions_; ///< Hold toroidal positions for toroidal scheme
    UnwrapType scheme_;   ///< Which unwrap scheme to use
};
#endif
