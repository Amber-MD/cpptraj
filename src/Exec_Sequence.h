#ifndef INC_EXEC_SEQUENCE_H
#define INC_EXEC_SEQUENCE_H
#include "Exec.h"
/// Create a molecule from 2 or more units 
class Exec_Sequence : public Exec {
  public:
    Exec_Sequence() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Sequence(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;
};
#endif
