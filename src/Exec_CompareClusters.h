#ifndef INC_EXEC_COMPARECLUSTERS_H
#define INC_EXEC_COMPARECLUSTERS_H
#include "Exec.h"
/// Compare two cluster number vs time sets 
class Exec_CompareClusters : public Exec {
  public:
    Exec_CompareClusters();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareClusters(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static DataSet* getClusterSet(std::string const&, DataSetList const&);
};
#endif
