#ifndef INC_EXEC_SET_H
#define INC_EXEC_SET_H
#include "Exec.h"
/// Used to set script variables
class Exec_Set : public Exec {
  public:
    Exec_Set() : Exec(CONTROL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Set(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static int AddVariable(CpptrajState&, std::string const&, std::string const&);
    static int AppendVariable(CpptrajState&, std::string const&, std::string const&);
    static int UpdateVariable(CpptrajState&, std::string const&, std::string const&);
};
#endif
