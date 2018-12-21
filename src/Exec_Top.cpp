#include "Exec_Top.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"
#include "TopInfo.h"

void Exec_LoadParm::Help() const {
  mprintf("\t<filename> [{[TAG] | name <setname>}]\n"
          "\t [{ nobondsearch |\n"
          "\t    [bondsearch <offset>] [searchtype {grid|pairlist}]\n"
          "\t  }]\n"
          "  Add <filename> to the topology list.\n"
          "  For topologies that may not have bond information, 'bondsearch <offset>'\n"
          "  controls the offset that will be added to atom-atom distances when\n"
          "  searching for bonds (default 0.2 Ang), and 'searchtype' specifies\n"
          "  different algorithms that can be used for searching bonds (still\n"
          "  experimental. Bond searching can be skipped via 'nobondsearch' (not\n"
          "  recommended).\n"
          "  Use 'help Formats parm' for format-specific options.\n");
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
static int CommonSetup(TopInfo& info, CpptrajState& State, ArgList& argIn, const char* desc)
{
  Topology* parm = 0;
  ReferenceFrame REF = State.DSL().GetReferenceFrame( argIn );
  if (REF.error()) return 1;
  if (REF.empty()) {
    parm = State.DSL().GetTopByIndex( argIn );
    if (parm == 0) return 1;
  } else
    mprintf("\tUsing '%s'\n", REF.refName());
  std::string outname = argIn.GetStringKey("out");
  int err = 0;
  if (outname.empty())
    err = info.SetupTopInfo( parm, REF.RefPtr() );
  else {
    CpptrajFile* outfile = State.DFL().AddCpptrajFile(outname, desc);
    if (outfile == 0) return CpptrajState::ERR;
    mprintf("\tOutput to '%s'\n", outfile->Filename().full());
    err = info.SetupTopInfo( outfile, parm, REF.RefPtr() );
  }
  return err;
}

// -----------------------------------------------------------------------------
void Exec_BondInfo::Help() const {
  mprintf("\t[%s] [<mask1>] [<mask2>] [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  For specified topology (first by default) either print bond info for all\n"
          "  atoms in <mask1>, or print info for bonds with first atom in <mask1> and\n"
          "  second atom in <mask2>.\n");
}

Exec::RetType Exec_BondInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Bond info")) return CpptrajState::ERR;
  std::string mask1 = argIn.GetMaskNext();
  if (info.PrintBondInfo( mask1, argIn.GetMaskNext(),false )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_UBInfo::Help() const {
  mprintf("\t[%s] [<mask1>] [<mask2>] [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  For specified topology (first by default) either print CHARMM Urey-Bradley\n"
          "  info for all atoms in <mask1>, or print info for bonds with first atom in\n"
          "  <mask1> and second atom in <mask2>.\n");
}

Exec::RetType Exec_UBInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Urey-Bradley info")) return CpptrajState::ERR;
  std::string mask1 = argIn.GetMaskNext();
  if (info.PrintBondInfo( mask1, argIn.GetMaskNext(), true )) return CpptrajState::ERR;
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
void Exec_AngleInfo::Help() const {
  mprintf("\t[%s] [<mask1>] [<mask2> <mask3>]\n\t[out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  For specified topology (first by default) either print angle info for all\n"
          "  atoms in <mask1>, or print info for angles with first atom in <mask1>,\n"
          "  second atom in <mask2>, and third atom in <mask3>.\n");
}

Exec::RetType Exec_AngleInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Angle info")) return CpptrajState::ERR;
  std::string mask1 = argIn.GetMaskNext();
  std::string mask2 = argIn.GetMaskNext();
  if (info.PrintAngleInfo( mask1, mask2, argIn.GetMaskNext() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_DihedralInfo::Help() const {
  mprintf("\t[%s] [<mask1>] [<mask2> <mask3> <mask4>]\n\t[out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  For specified topology (first by default) either print dihedral info for all\n"
          "  atoms in <mask1>, or print info for dihedrals with first atom in <mask1>,\n"
          "  second atom in <mask2>, third atom in <mask3>, and fourth atom in <mask4>.\n");
}

Exec::RetType Exec_DihedralInfo::Execute(CpptrajState& State, ArgList& argIn) {
  if (argIn.hasKey("and")) {
    mprinterr("Error: The 'and' keyword has been deprecated. To restrict dihedral\n"
              "Error:   selection please use 4 masks.\n");
    return CpptrajState::ERR;
  }
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Dihedral info")) return CpptrajState::ERR;
  std::string mask1 = argIn.GetMaskNext();
  std::string mask2 = argIn.GetMaskNext();
  std::string mask3 = argIn.GetMaskNext();
  if (info.PrintDihedralInfo( mask1, mask2, mask3, argIn.GetMaskNext(), false ))
    return CpptrajState::ERR;
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
void Exec_ImproperInfo::Help() const {
  mprintf("\t[%s] [<mask1>] [<mask2> <mask3> <mask4>]\n\t[out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  For specified topology (first by default) either print CHARMM improper info\n"
          "  for all atoms in <mask1>, or print info for dihedrals with first atom in <mask1>,\n"
          "  second atom in <mask2>, third atom in <mask3>, and fourth atom in <mask4>.\n");
}

Exec::RetType Exec_ImproperInfo::Execute(CpptrajState& State, ArgList& argIn) {
  if (argIn.hasKey("and")) {
    mprinterr("Error: The 'and' keyword has been deprecated. To restrict improper\n"
              "Error:   selection please use 4 masks.\n");
    return CpptrajState::ERR;
  }
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Improper info")) return CpptrajState::ERR;
  std::string mask1 = argIn.GetMaskNext();
  std::string mask2 = argIn.GetMaskNext();
  std::string mask3 = argIn.GetMaskNext();
  if (info.PrintDihedralInfo( mask1, mask2, mask3, argIn.GetMaskNext(), true ))
    return CpptrajState::ERR;
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
void Exec_AtomInfo::Help() const {
  mprintf("\t[%s] [<mask>] [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print information on atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_AtomInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Atom info")) return CpptrajState::ERR;
  if (info.PrintAtomInfo( argIn.GetMaskNext() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}

// -----------------------------------------------------------------------------
void Exec_ResInfo::Help() const {
  mprintf("\t[%s] [<mask>] [short [maxwidth <#res>]]\n\t[out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for residues in <mask> for specified topology (first by default).\n"
          "  If 'short' is specified print residue info in shorter form.\n");
}

Exec::RetType Exec_ResInfo::Execute(CpptrajState& State, ArgList& argIn) {
  bool printShort = argIn.hasKey("short");
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Residue info")) return CpptrajState::ERR;
  int err;
  if (printShort)
    err = info.PrintShortResInfo( argIn.GetMaskNext(), argIn.getKeyInt("maxwidth",50) );
  else
    err = info.PrintResidueInfo( argIn.GetMaskNext() );
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_MolInfo::Help() const {
  mprintf("\t[%s] [<mask>] [short] [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print info for molecules in <mask> for specfied topology (first by default).\n"
          "  If 'short' is specified print counts of each unique molecule.\n");
}

Exec::RetType Exec_MolInfo::Execute(CpptrajState& State, ArgList& argIn) {
  bool printShort = argIn.hasKey("short");
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Molecule info")) return CpptrajState::ERR;
  int err;
  if (printShort)
    err = info.PrintShortMolInfo( argIn.GetMaskNext() );
  else
    err = info.PrintMoleculeInfo( argIn.GetMaskNext() );
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_ChargeInfo::Help() const {
  mprintf("\t[%s] <mask> [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print total charge of atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_ChargeInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Charge info")) return CpptrajState::ERR;
  if (info.PrintChargeInfo( argIn.GetMaskNext() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_MassInfo::Help() const {
  mprintf("\t[%s] <mask> [out <file>]\n", DataSetList::TopIdxArgs);
  mprintf("  Print total mass of atoms in <mask> for specified topology (first by default).\n");
}

Exec::RetType Exec_MassInfo::Execute(CpptrajState& State, ArgList& argIn) {
  TopInfo info;
  if (CommonSetup(info, State, argIn, "Mass info")) return CpptrajState::ERR;
  if (info.PrintMassInfo( argIn.GetMaskNext() )) return CpptrajState::ERR;
  return CpptrajState::OK;
}
