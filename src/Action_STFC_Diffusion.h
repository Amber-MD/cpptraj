#ifndef INC_ACTION_STFC_DIFFUSION_H
#define INC_ACTION_STFC_DIFFUSION_H
#include "Action.h"
#include "ImagedAction.h"
/** \author Hannes H. Loeffler
  * \author C++ adaptation by Daniel R. Roe
  */
class Action_STFC_Diffusion : public Action {
  public:
    Action_STFC_Diffusion();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_STFC_Diffusion(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    void calculateMSD(const double*,int,int,Vec3 const&);

    bool printDistances_; // iarg1
    enum CalcType { DEFAULT = 0, COM, DIST };
    CalcType calcType_; // iarg2
    enum DirectionType { DX = 0, DY, DZ, DXY, DXZ, DYZ, DXYZ };
    DirectionType direction_; // iarg3
    AtomMask mask_;
    AtomMask mask2_;
    CpptrajFile* output_;
    CpptrajFile* outputnw_;
    CpptrajFile* outputad_;
    double time_;
    double lowerCutoff_;
    double upperCutoff_;
    bool hasBox_;
    int n_atom_;

    typedef std::vector<double> Darray;
    Darray initialxyz_;
    Darray distancexyz_;
    Darray distance_;
    Darray deltaxyz_;
    Darray previousxyz_;

    Darray dSum1_;
    Darray dSum2_;
    std::vector<int> nInside_;
    int elapsedFrames_;
    ImagedAction image_; ///< Imaging routines.
};
#endif    
