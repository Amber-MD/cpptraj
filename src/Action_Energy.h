#ifndef INC_ACTION_ENERGY_H
#define INC_ACTION_ENERGY_H
#include "Action.h"
#include "Energy.h"
/// Calculate energy 
class Action_Energy: public Action {
  public:
    Action_Energy();
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
    /// Corresponds to calculations.
    enum CalcType { BND, ANG, DIH, N14, NBD, LJ, COULOMB };
    std::vector<DataSet*> Energy_; ///< Hold output data sets
    std::vector<CalcType> Ecalcs_; ///< Hold which calcs to perform
    typedef std::vector<CalcType>::const_iterator calc_it;
    Topology* currentParm_;        ///< Hold current topology
    CharMask Mask1_;               ///< Char mask for all but NB calc
    AtomMask Imask_;               ///< Int mask for NB calc
    Energy_Amber ENE_;             ///< Energy calc class.
};
#endif
