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

    double E_bond(Frame const&, Topology const&, CharMask const&);
    double E_angle(Frame const&, Topology const&, CharMask const&);
    double E_torsion(Frame const&, Topology const&, CharMask const&);
    double E_14_Nonbond(Frame const&, Topology const&, CharMask const&, double&);
    double E_Nonbond(Frame const&, Topology const&, AtomMask const&, double&);

    void SetDebug(int d) { debug_ = d; }
    void PrintTiming() const;
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
    Timer time_bond_;
    Timer time_angle_;
    Timer time_tors_;
    Timer time_14_;
    Timer time_NB_;
};
#ifdef USE_SANDERLIB
class Energy_Sander {
  public:
    Energy_Sander() : top_pindex_(-1) {}
    ~Energy_Sander();
    int Initialize(Topology const&, Frame&);
    int CalcEnergy(Topology const&, Frame&);
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
    FileName top_filename_;
    std::vector<double> forces_;
    int top_pindex_;
};
#endif /* USE_SANDERLIB */
#endif 
