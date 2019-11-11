#ifndef INC_ACTION_ENERGY_H
#define INC_ACTION_ENERGY_H
#include "Action.h"
#include "Energy.h"
#include "CharMask.h"
#include "Timer.h"
class Ewald;
/// Calculate energy 
class Action_Energy: public Action {
  public:
    Action_Energy();
    ~Action_Energy();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Energy(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
    /// Corresponds to data sets.
    enum Etype { BOND = 0, ANGLE, DIHEDRAL, V14, Q14, VDW, ELEC, KE, TOTAL};
    /// Add energy data set of specified type.
    int AddSet(Etype, DataSetList&, DataFile*, std::string const&);
    /// For debugging the direct sum convergence
    double Dbg_Direct(Frame const&,int);
    /// Corresponds to calculations.
    enum CalcType { C_BND = 0, C_ANG, C_DIH, C_N14, C_NBD, C_LJ,
                    C_COULOMB, C_DIRECT, C_EWALD, C_PME, C_KEAUTO, C_KEVEL, C_KEVV };
    /// Corresponds to type of electrostatics.
    enum ElecType { NO_ELE = 0, SIMPLE, DIRECTSUM, EWALD, PME };
    /// Corresponds to type of KE calc.
    enum KEType { KE_NONE = 0, KE_AUTO, KE_VEL, KE_VV };

    ElecType elecType_;            ///< Type of electrostatics calc.
    KEType KEtype_;                ///< Type of KE calc.
    std::vector<DataSet*> Energy_; ///< Hold output data sets (length Etype+1)
    std::vector<CalcType> Ecalcs_; ///< Hold which calcs to perform
    Topology* currentParm_;        ///< Hold current topology
    CharMask Mask1_;               ///< Char mask for all but NB calc
    AtomMask Imask_;               ///< Int mask for NB calc
    Energy_Amber ENE_;             ///< Energy calc class.
    std::string setname_;          ///< Output DataSet name
    int npoints_;                  ///< # cells in each direction (DIRECT) or spline order (PME)
    int debug_;
    Ewald* EW_;                    ///< Ewald energy class.
    double cutoff_;                ///< Ewald direct space cutoff.
    double dsumtol_;               ///< Ewald direct sum tolerance.
    double rsumtol_;               ///< Regular Ewald reciprocal sum tolerance.
    double ewcoeff_;               ///< Ewald coefficient.
    double lwcoeff_;               ///< LJ Ewald coefficient.
    double ljswidth_;              ///< Size of LJ switch region
    double maxexp_;
    double skinnb_;                ///< Size of non-bonded "skin"
    double erfcDx_;                ///< Spacing for ERFC table (default 1/5000)
    double dt_;                    ///< Time step for estimating kinetic energy (leapfrog)
    int mlimits_[3];               ///< mlimits (reg. Ewald) or nfft (PME)
    bool need_lj_params_;          ///< True if LJ parameters needed.
    Timer time_total_;
    Timer time_bond_;
    Timer time_angle_;
    Timer time_tors_;
    Timer time_14_;
    Timer time_NB_;
    Timer time_ke_;
};
#endif
