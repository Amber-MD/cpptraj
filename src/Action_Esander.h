#ifndef INC_ACTION_ESANDER_H
#define INC_ACTION_ESANDER_H
#include "Action.h"
#ifdef USE_SANDERLIB
#  include "Energy_Sander.h"
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
    /// Add energy data set of specified type.
    inline int AddSet(Energy_Sander::Etype, DataSetList&, DataFile*, std::string const&);
    /// Set up sander energy
    int InitForRef();
    typedef std::vector<DataSet*> Earray;
    Earray Esets_;          ///< Hold output data sets
    Frame refFrame_;        ///< Hold reference coords for init.
    Topology* currentParm_; ///< Hold current topology
    Energy_Sander SANDER_;  ///< Sander energy routines.
    bool save_forces_;      ///< If true save force information.
    CoordinateInfo cInfo_;  ///< Coordinate info if saving force info.
    Frame newFrame_;        ///< Frame for saving force info.
    Action::RetType ret_;   ///< Return status for DoAction()
    ActionInit Init_;       ///< Hold master DataSet/DataFile lists.
    std::string setname_;   ///< Data set name.
    DataFile* outfile_;     ///< Output data file.
#   endif /* USE_SANDERLIB */
};
#endif
