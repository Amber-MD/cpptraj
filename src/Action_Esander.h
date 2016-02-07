#ifndef INC_ACTION_ESANDER_H
#define INC_ACTION_ESANDER_H
#include "Action.h"
#ifdef USE_SANDERLIB
#  include "Esander.h"
#endif
/// Calculate energy via sanderlib
class Action_Esander: public Action {
  public:
    Action_Esander();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Esander(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef USE_SANDERLIB
    /// Corresponds to data sets.
    enum Etype { BOND = 0, ANGLE, DIHEDRAL, V14, Q14, VDW, ELEC, TOTAL};
    /// Add energy data set of specified type.
    int AddSet(Etype, DataSetList&, DataFile*, std::string const&);
    typedef std::vector<DataSet*> Earray;
    Earray Esets_;          ///< Hold output data sets
    Frame refFrame_;        ///< Hold reference coords for init.
    Topology* currentParm_; ///< Hold current topology
    Energy_Sander SANDER_;  ///< Sander energy routines.
#   endif /* USE_SANDERLIB */
};
#endif
