#ifndef INC_ACTION_ENERGY_H
#define INC_ACTION_ENERGY_H
#include "Action.h"
#include "Energy.h"
#include "Ewald.h"
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
    enum Etype { BOND = 0, ANGLE, DIHEDRAL, V14, Q14, VDW, ELEC, TOTAL};
    /// Add energy data set of specified type.
    int AddSet(Etype, DataSetList&, DataFile*, std::string const&);
    /// For debugging the direct sum convergence
    double Dbg_Direct(Frame const&,int);
    /// Corresponds to calculations.
    enum CalcType { C_BND, C_ANG, C_DIH, C_N14, C_NBD, C_LJ,
                    C_COULOMB, C_DIRECT, C_EWALD, C_PME };
    /// Corresponds to type of electrostatics.
    enum ElecType { NO_ELE = 0, SIMPLE, DIRECTSUM, EWALD, PME };

    ElecType elecType_;            ///< Type of electrostatics calc.
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
    double maxexp_;
    double skinnb_;                ///< Size of non-bonded "skin"
    double erfcDx_;                ///< Spacing for ERFC table (default 1/5000)
    int mlimits_[3];               ///< mlimits (reg. Ewald) or nfft (PME)
    bool need_lj_params_;          ///< True if LJ parameters needed.
    Timer etime_;
};
#endif
