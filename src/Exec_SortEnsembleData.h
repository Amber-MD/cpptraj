#ifndef INC_EXEC_SORTENSEMBLEDATA_H
#define INC_EXEC_SORTENSEMBLEDATA_H
#include "Exec.h"
/// <Enter description of Exec_SortEnsembleData here>
class Exec_SortEnsembleData : public Exec {
  public:
    Exec_SortEnsembleData() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SortEnsembleData(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int SortData(DataSetList const&) const;
    int Sort_pH_Data(DataSetList const&) const;
};
#endif
