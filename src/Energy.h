#ifndef INC_ENERGY_H
#define INC_ENERGY_H
#include "ParameterTypes.h"
class Topology;
class AtomMask;
class CharMask;
class Frame;
class ExclusionArray;
/// Calculate energy/force from coordinates.
class Energy_Amber {
  public:
    typedef std::vector<double> Darray;
    Energy_Amber();

    double E_bond(Frame const&, Topology const&, CharMask const&);
    double E_angle(Frame const&, Topology const&, CharMask const&);
    double E_torsion(Frame const&, Topology const&, CharMask const&);
    double E_14_Nonbond(Frame const&, Topology const&, CharMask const&, double&);
    double E_Nonbond(Frame const&, Topology const&, AtomMask const&, double&, ExclusionArray const&);
    double E_VDW(Frame const&, Topology const&, AtomMask const&, ExclusionArray const&);
    double E_Elec(Frame const&, Topology const&, AtomMask const&, ExclusionArray const&);

    double E_DirectSum(Frame const&, Topology const&, AtomMask const&, ExclusionArray const&,int);
    /// Calculate kinetic energy from velocity information.
    double E_Kinetic(Frame const&, AtomMask const&);
    /// Calculate kinetic energy from forces and plus-half timestep velocities.
    double E_Kinetic_VV(Frame const&, AtomMask const&, double);

    void SetDebug(int d) { debug_ = d; }
  private:
    double CalcBondEnergy(Frame const&, BondArray const&, BondParmArray const&,
                          CharMask const&);
    double CalcAngleEnergy(Frame const&, AngleArray const&, AngleParmArray const&,
                           CharMask const&);
    double CalcTorsionEnergy(Frame const&, DihedralArray const&, DihedralParmArray const&,
                             CharMask const&);
    double Calc_14_Energy(Frame const&, DihedralArray const&, DihedralParmArray const&,
                          Topology const&, CharMask const&, double&);

    static const double QFAC;
    int debug_;
};
#endif 
