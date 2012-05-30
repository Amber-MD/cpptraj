#ifndef INC_ACTION_CONTACTS_H
#define INC_ACTION_CONTACTS_H
#include <map>
#include <set>
#include "Action.h"
class Action_Contacts : public Action {
  public:
    Action_Contacts();
    ~Action_Contacts();
  private:
    AtomMask Mask_;
    bool byResidue_;
    double distance_;
    //Topology* RefParm_;
    bool first_;
    CpptrajFile outfile_;
    CpptrajFile outfile2_;
    typedef std::pair<int,int> contactType;
    //typedef std::vector< contactType > contactListType;
    typedef std::multimap<int,int> contactListType;
    contactListType nativecontacts_;
    std::vector<int> residueContacts_;
    std::vector<int> residueNative_;
    std::set<int> activeResidues_;

    int init();
    int setup();
    int action();

    int SetupContacts(Frame*, Topology*);
};
#endif
