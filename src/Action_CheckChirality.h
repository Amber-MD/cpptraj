#ifndef INC_ACTION_CHECKCHIRALITY_H
#define INC_ACTION_CHECKCHIRALITY_H
#include "Action.h"
#include "Array1D.h"
/// Determine if amino acids are D or L 
class Action_CheckChirality: public Action {
  public:
    Action_CheckChirality() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_CheckChirality(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
#   endif
    void Print();

    struct ResidueInfo {
//      DataSet* data_; ///< data
      int num_; ///< residue number
      int isActive_;
      int n_;   ///< N coord index
      int ca_;  ///< C alpha coord index
      int c_;   ///< C coord index
      int cb_;  ///< C beta coord index
      int N_L_; ///< # times this residue was L
      int N_D_; ///< # times this residue was D
    };

    typedef std::vector<ResidueInfo> Rarray;
    Rarray resInfo_; 
    CharMask Mask1_;
    DataSet* data_L_;     ///< Hold number of times each residue was L
    DataSet* data_D_;     ///< Hold number of times each residue was D
    std::string setname_; ///< Data set name
    ActionInit Init_;     ///< Master DSL/DFL
};
#endif
