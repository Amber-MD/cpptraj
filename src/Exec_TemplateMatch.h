#ifndef INC_EXEC_TEMPLATEMATCH_H
#define INC_EXEC_TEMPLATEMATCH_H
#include "Exec.h"
/// Immediate command: map a target topology onto a user template for TI.
/** Alias: timap. Source encoding is UTF-8 (comments may use O5′, χ, λ).
  * Looks up <tgt> / template as a TOPOLOGY set, then a COORDS set (Amber
  * OFF units loaded with readdata appear as Name[Unit], e.g. FLE[FLE]),
  * then a numeric parm index. See TemplateMatch for the matching algorithm.
  */
class Exec_TemplateMatch : public Exec {
  public:
    Exec_TemplateMatch() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_TemplateMatch(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
