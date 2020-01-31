#ifndef INC_EXEC_HELP_H
#define INC_EXEC_HELP_H
#include "Exec.h"
/// Find help for command/topic
class Exec_Help : public Exec {
  public:
    Exec_Help() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Help(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int Formats(ArgList&) const;
    int Masks(ArgList&) const;
    int Math(ArgList&) const;
    int Topics(ArgList&) const;
};
#endif
