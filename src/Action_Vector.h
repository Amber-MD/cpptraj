#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "Structure/LeastSquaresPlane.h" // For CORRPLANE
class DataSet_Vector;
class DataSet_3D;
class CharMask;
class Action_Vector : public Action {
  public:
    Action_Vector();
    ~Action_Vector();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Vector(); }
    void Help() const;
  private:
    enum vectorMode {
      NO_OP=0,   PRINCIPAL_X, PRINCIPAL_Y, PRINCIPAL_Z,
      DIPOLE,    BOX,         MASK,
      CORRPLANE, CENTER,      BOX_X,       BOX_Y,       BOX_Z,
      BOX_CTR,   MINIMAGE,    MOMENTUM,    VELOCITY,    FORCE,
      BONDDIPOLE
    };
    static const char* ModeString_[];
    static const bool NeedsOrigin_[];

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    /// \return Center of mass or geometric center of atoms in given mask
    inline Vec3 GetVec(Frame const&, AtomMask const&) const;
    void Mask(Frame const&);
    static inline bool calcBondDipole(Vec3&, Vec3&, Vec3&, double, Vec3 const&, Vec3 const&);
    void BondDipole_individualBonds(Frame const&);
    void BondDipole_net_bondOrigin(Frame const&);
    void Dipole(Frame const&);
    void Principal(Frame const&);
    void CorrPlane(Frame const&);
    void UnitCell(Box const&, Vec3 const&);
    void BoxLengths(Box const&);
    void MinImage(Frame const&);

    DataSet_Vector* Vec_;   ///< Hold vector values
    DataSet* Magnitude_;    ///< Hold vector magnitudes if requested
    DataSet_3D* gridSet_;   ///< Hold grid set for getting box vectors from grid.
    Cpptraj::Structure::LeastSquaresPlane vcorr_; ///< Temp. space for calculating CorrPlane
    vectorMode mode_;       ///< Vector calculation mode
    bool ptrajoutput_;      ///< If true output in ptraj format
    bool needBoxInfo_;      ///< If true box info required.
    bool useMass_;          ///< If true, centers are mass-weighted 
    bool dipole_in_debye_;  ///< If true, report dipole vector values in Debye
    Topology* CurrentParm_; ///< Current topology (for dipole)
    AtomMask mask_;
    AtomMask mask2_;
    CpptrajFile* outfile_;
    CharMask* cmask_;
    int debug_;
};
#endif
