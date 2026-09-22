#ifndef INC_EXEC_TIMAP_H
#define INC_EXEC_TIMAP_H
#include "Exec.h"
/// Immediate command: write a residue library in template atom order; optional TI pair.
/** Command: timap (alias: timatch).
  * Resolves <tgt> / template / series members to Topology via:
  *   1. Named TOPOLOGY set (parm mol2/pdb/amber/...),
  *   2. Named COORDS set (Amber OFF .lib via readdata → Name[Unit]),
  *   3. Numeric parm index,
  *   4. Or a file path — auto-loaded with ParmFile (mol2, pdb, ...) or
  *      DataFile Amber OFF reader (.lib), then matched as Topology.
  * Formats may be mixed in one call; TIMatch always sees Topology.
  * Default: write the remapped target (Amber OFF .lib, mol2, or pdb via
  * out / fmt). Optional tiout writes dual-topology mol2/lib (same NATOM,
  * charge-0 / mass-0 / type DUM dummies) plus scmask for pmemd.
  * naorder / aaorder walk <tgt> itself when no template is given.
  * series mode auto-selects a parent among many analogs (leadopt-inspired).
  * See TIMatch for the matching algorithm.
  * \author Nathan D. Levinzon <ndlevinzon@gmail.com>
  */
class Exec_TIMap : public Exec {
  public:
    Exec_TIMap() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_TIMap(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
