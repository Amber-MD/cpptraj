#ifndef INC_ENERGY_H
#define INC_ENERGY_H
#include "Topology.h"
#include "Timer.h"
#ifdef USE_SANDERLIB
#  include "sander.h"
#endif
/// Calculate energy/force from coordinates.
class Energy_Amber {
  public:
    typedef std::vector<double> Darray;
    Energy_Amber();

    double E_bond(Frame const&, Topology const&, AtomMask const&);
    double E_angle(Frame const&, Topology const&, AtomMask const&);
    double E_torsion(Frame const&, Topology const&, AtomMask const&);
    double E_14_Nonbond(Frame const&, Topology const&, AtomMask const&, double&);
    double E_Nonbond(Frame const&, Topology const&, AtomMask const&, double&);

    void SetDebug(int d) { debug_ = d; }
    void PrintTiming() const;
  private:
    double CalcBondEnergy(Frame const&, BondArray const&, BondParmArray const&,
                          AtomMask const&);
    double CalcAngleEnergy(Frame const&, AngleArray const&, AngleParmArray const&,
                           AtomMask const&);
    double CalcTorsionEnergy(Frame const&, DihedralArray const&, DihedralParmArray const&,
                             AtomMask const&);
    double Calc_14_Energy(Frame const&, DihedralArray const&, DihedralParmArray const&,
                          Topology const&, AtomMask const&, double&);

    static const double QFAC;
    int debug_;
    Timer time_bond_;
    Timer time_angle_;
    Timer time_tors_;
    Timer time_14_;
    Timer time_NB_;
};
class Energy_Sander {
  public:
#ifdef USE_SANDERLIB
    Energy_Sander() : top_(0) {}
    ~Energy_Sander();
    int Initialize(Topology*, Frame&);
    int CalcEnergy(Topology*, Frame&);
    double Ebond() const { return energy_.bond; }
    double Eangle() const { return energy_.angle; }
    double Edihedral() const { return energy_.dihedral; }
    double Evdw14() const { return energy_.vdw_14; }
    double Eelec14() const { return energy_.elec_14; }
    double Evdw() const { return energy_.vdw; }
    double Eelec() const { return energy_.elec; } 
  private:
    sander_input input_;
    pot_ene energy_;
    Topology* top_;
    std::vector<double> forces_;
#else
    Energy_Sander() {}
    int Initialize(Topology*,Frame&) { return 1; }
    int CalcEnergy(Topology*,Frame&) { return 1; }
    double Ebond() const { return 0.0; }
    double Eangle() const { return 0.0; }
    double Edihedral() const { return 0.0; }
    double Evdw14() const { return 0.0; }
    double Eelec14() const { return 0.0; }
    double Evdw() const { return 0.0; }
    double Eelec() const { return 0.0; }
#endif 
};
#endif 
