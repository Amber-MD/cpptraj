#ifndef INC_EXEC_PARSETIMING_H
#define INC_EXEC_PARSETIMING_H
#include "Exec.h"
/// For extracting and parsing timing data from cpptraj output 
class Exec_ParseTiming : public Exec {
  public:
    Exec_ParseTiming();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParseTiming(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    class RunTiming;

    RunTiming read_cpptraj_output(std::string const&);
};
#endif
