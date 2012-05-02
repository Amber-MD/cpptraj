#ifndef INC_ACTION_CONTACTS_H
#define INC_ACTION_CONTACTS_H
#include "Action.h"
//#include "ActionReference.h"
// TODO: Use with ActionReference?
class Action_Contacts : public Action {
  public:
    Action_Contacts();
    ~Action_Contacts();
  private:
    AtomMask Mask_;
    bool byResidue_;
    double distance_;
    Topology* RefParm_;
    bool first_;
    CpptrajFile outfile_;
    std::vector< std::pair<int,int> > contactlist_;

    int init();
    int setup();
    int action();

    int SetupContacts(Frame*);
};
#endif
