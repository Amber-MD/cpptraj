#ifndef INC_ACTION_CONTACTS_H
#define INC_ACTION_CONTACTS_H
#include <map>
#include <set>
#include "Action.h"
class Action_Contacts : public Action {
  public:
    Action_Contacts();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Contacts(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    int SetupContacts(Frame const&, Topology const&);

    AtomMask Mask_;
    bool byResidue_;
    double distance_;
    double dt_;
    bool first_;
    Topology* CurrentParm_;
    CpptrajFile* outfile_;
    CpptrajFile* outfile2_;
    typedef std::pair<int,int> contactType;
    //typedef std::vector< contactType > contactListType;
    typedef std::multimap<int,int> contactListType;
    contactListType nativecontacts_;
    std::vector<int> residueContacts_;
    std::vector<int> residueNative_;
    std::set<int> activeResidues_;
};
#endif
