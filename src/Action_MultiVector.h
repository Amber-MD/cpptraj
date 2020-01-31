#ifndef INC_MULTIVECTOR_H
#define INC_MULTIVECTOR_H
#include "Action.h"
#include "Range.h"
class DataSet_Vector;
class Action_MultiVector : public Action {
  public:
    Action_MultiVector();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MultiVector(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    int debug_;
    std::vector<DataSet_Vector*> data_;  ///< Output DataSets, 1 per vector.
    Range resRange_;                     ///< Residues to search for vectors.
    std::string dsetname_;
    NameType name1_;
    NameType name2_;
    std::vector<int> CrdIdx1_;
    std::vector<int> CrdIdx2_;
    DataFile* outfile_;
    // TODO: Replace these with new DataSet type
    DataSetList* masterDSL_;
    bool ired_;
};
#endif 
