#ifndef INC_EXEC_ANALYZE_H
#define INC_EXEC_ANALYZE_H
#include "Exec.h"
/// Add an analysis command to the analysis queue. For backwards compat.
class Exec_Analyze : public Exec {
  public:
    Exec_Analyze() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Analyze(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
