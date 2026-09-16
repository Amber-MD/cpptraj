#ifndef INC_EXEC_TIMATCH_H
#define INC_EXEC_TIMATCH_H
#include "Exec.h"
/// Immediate command: write a residue library in template atom order; optional TI pair.
/** Command: timatch (aliases: templatematch, timap).
  * Looks up <tgt> / template as a TOPOLOGY set, then a COORDS set (Amber
  * OFF units loaded with readdata appear as Name[Unit], e.g. FLE[FLE]),
  * then a numeric parm index.
  * Default: write the target as an Amber OFF .lib in template atom order
  * (out <file>; default <residue>.sorted.lib). No dummy atoms.
  * Optional tiout <prefix> writes dual-topology mol2/lib (same NATOM,
  * charge-0 / mass-0 / type DUM dummies) plus scmask for pmemd.
  * See TemplateMatch for the matching algorithm.
  */
class Exec_TIMatch : public Exec {
  public:
    Exec_TIMatch() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_TIMatch(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
