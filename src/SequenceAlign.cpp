#include "SequenceAlign.h"
// EXPERIMENTAL ALPHA CODE
/** Idea is to take a reference PDB and a BLAST-type alignment
  * with format:
  * Query  1    MIVNVIQKDRL-KEQKLQFIRNHQQAFDVEPIRPLPLFEDFVTSIEGDCSLEASCKIESD  59
  *             MI+     D   KE+++Q +R+H ++FDVE   PLPLFE  V S++    LE S K++
  * Sbjct  1    MIMTTTWPDSYAKERRIQRLRHHFESFDVERAFPLPLFEQAVLSLDSCPLLEPSFKVQEG  60
  * and split the reference PDB into segments with the correct residue IDs so
  * it can then be put together.
  *
  * Usage: sequencealign <ref keyword> blastfile <file>
  */
int SequenceAlign(CpptrajState& State, ArgList& argIn) {

  return 0;
}
