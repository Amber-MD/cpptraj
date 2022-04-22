#ifndef INC_EXEC_COMPAREENERGY_H
#define INC_EXEC_COMPAREENERGY_H
#include "Exec.h"
/// <Enter description of Exec_CompareEnergy here>
class Exec_CompareEnergy : public Exec {
  public:
    Exec_CompareEnergy() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareEnergy(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static DataSet_Coords* GetCoordsSet(DataSetList const&, std::string const&);
};
#endif
