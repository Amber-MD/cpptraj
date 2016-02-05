#include "Exec_Top.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

void Exec_LoadParm::Help() const {
  mprintf("\t<filename> [{[TAG] | name <setname>}] [nobondsearch | bondsearch [<offset>]]\n"
          "  Add <filename> to the topology list.\n");
  ParmFile::ReadOptions();
}
// -----------------------------------------------------------------------------
void Exec_ParmInfo::Help() const {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print information on specfied topology (first by default).\n");
}

Exec::RetType Exec_ParmInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->Summary();
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_BondInfo::Help() const {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print bond info for atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_BondInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->PrintBondInfo( argIn.GetMaskNext() );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_AngleInfo::Help() const {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print angle info for atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_AngleInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->PrintAngleInfo( argIn.GetMaskNext() );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_DihedralInfo::Help() const {
  mprintf("\t[%s] [<mask>] [and]\n", DataSetList::TopIdxArgs);
  mprintf("  Print dihedral info for atoms in <mask> for specified topology (first by default).\n"
          "  If 'and' is specified dihedral must include all atoms specfied by <mask>.\n");
}

Exec::RetType Exec_DihedralInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->PrintDihedralInfo( argIn.GetMaskNext(), !argIn.hasKey("and") );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_AtomInfo::Help() const {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print information on atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_AtomInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->PrintAtomInfo( argIn.GetMaskNext() );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_ResInfo::Help() const {
  mprintf("\t[%s] [<mask>] [short]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for residues in <mask> for specified topology (first by default).\n"
          "  If 'short' is specified print residue info in shorter form.\n");
}

Exec::RetType Exec_ResInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  bool printShort = argIn.hasKey("short");
  if (printShort)
    parm->PrintShortResInfo( argIn.GetMaskNext(), argIn.getKeyInt("maxwidth",50) );
  else
    parm->PrintResidueInfo( argIn.GetMaskNext() );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_MolInfo::Help() const {
  mprintf("\t[%s] [<mask>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for molecules in <mask> for specfied topology (first by default).\n");
}

Exec::RetType Exec_MolInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->PrintMoleculeInfo( argIn.GetMaskNext() );
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_ChargeInfo::Help() const {
  mprintf("\t[%s] <mask>\n", DataSetList::TopIdxArgs);
  mprintf("  Print total charge of atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_ChargeInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 0 )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_MassInfo::Help() const {
  mprintf("\t[%s] <mask>\n", DataSetList::TopIdxArgs);
  mprintf("  Print total mass of atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_MassInfo::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 1 )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
