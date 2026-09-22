#include "Exec_TIMap.h"
#include "TIMatch.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "CpptrajState.h"
#include "DataSet_Topology.h"
#include "DataSet_Coords.h"
#include "DataFile.h"
#include "ParmFile.h"
#include "FileName.h"
#include "StringRoutines.h"
#include "ArgList.h"
#include "Frame.h"
#include "CoordinateInfo.h"
#include "Trajout_Single.h"
#include "Trajin_Single.h"
#include "TrajectoryFile.h"
#include "Residue.h"
#include "Atom.h"
#include "Topology.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <vector>

/* Source encoding: UTF-8. Unicode (O5′, 2′, λ) is comments-only.
 *
 * Workflow:
 *   1. Match target onto template (partial map).
 *   2. Write the remapped target (out / fmt: Amber OFF .lib, mol2, or pdb).
 *   3. If tiout <prefix> is given, also build a dual-topology pair of the
 *      same size: unique atoms are real in one end state and charge-0 /
 *      mass-0 / type-DUM dummies in the other.
 *   4. series mode: choose one parent among many analogs (leadopt-inspired
 *      pairwise score), then map every other analog onto that parent.
 *
 * Series parent selection cites ideas from util/leadopt/ (not linked into the
 * build). Authors / sources:
 *   - Jonathan Redmann & Christopher Summa, Summa Lab, University of New
 *     Orleans — GraphGenerator4.py (TI graph planning from similarity).
 *   - leadopt/similarity.py — exp(-BETA * delta) similarity from atom-count
 *     differences relative to a common substructure (BETA = 0.1).
 *   - leadopt/mcs.py, rule.py, graph.py — pairwise MCS / common-substructure
 *     scoring used to decide which compounds share the most.
 * Here we approximate the common substructure with TIMatch's partial
 * map (nMapped / insertions / unmatched) rather than an external MCS engine.
 * \author Nathan D. Levinzon <ndlevinzon@gmail.com>
 */

// Exec_TIMap::Help()
void Exec_TIMap::Help() const {
  mprintf("\t{ <tgt> [template <name>] [naorder|aaorder] |\n"
          "\t  series <a> <b> ... [outprefix <pfx>] [tioutprefix <pfx>]\n"
          "\t         [mapoutprefix <pfx>] [parentout <file>] }\n"
          "\t[out <file>] [fmt {lib|mol2|pdb}] [tiout <prefix>] [mapout <file>] [maponly]\n"
          "\t[name <newparm>] [replace]\n"
          "\t[seed {auto|names|na|aa|none}] [anchor <atomname>]\n"
          "  Reorder atoms so shared sites share indices for TI (nucleotides,\n"
          "  amino acids, small molecules). Partial maps are expected (unlike atommap).\n"
          "\n"
          "  Atom order (pick one; combine template with naorder/aaorder to walk the\n"
          "  parent first when its file order is messy):\n"
          "    template <name>  parent's current order (Amber .lib, mol2, ...).\n"
          "    naorder          Amber NA walk of <tgt> (or of template if given):\n"
          "                     P -> OP -> O5' -> C5' -> C4' -> O4' -> C1' -> base ->\n"
          "                     C3' -> C2' -> O3'.\n"
          "    aaorder          ff19SB amino19.lib walk of <tgt> (or of template):\n"
          "                     N -> H -> CA -> HA -> side chain from CB -> C -> O\n"
          "                     (proline: N -> CD -> ... -> CB -> CA -> C -> O).\n"
          "    series <a> <b>... auto-pick one parent among the listed topologies\n"
          "                     (max mapped-atom sum; ties: fewest atoms, then\n"
          "                     leadopt-style exp(-0.1*(ins+unmap)) sum; see manual),\n"
          "                     then map every member onto it. Not with template.\n"
          "\n"
          "  Output:\n"
          "    Default: remapped structure (out <file>; else <residue>.sorted.lib).\n"
          "    fmt {lib|mol2|pdb}  choose writer (default lib). If out has extension\n"
          "                       .mol2 / .pdb / .lib, that format is used instead.\n"
          "    maponly  skip the aligned structure / in-memory parm; still write\n"
          "             mapout/tiout.\n"
          "    mapout   atom correspondence file.\n"
          "    name / replace  store remapped topology, or overwrite <tgt>.\n"
          "    tiout <prefix>  dual-topology TI pair (same NATOM; DUM q=0 m=0):\n"
          "      <prefix>.0.mol2/.lib  lambda=0 template + dummy insertions\n"
          "      <prefix>.1.mol2/.lib  lambda=1 target + dummy unmatched atoms\n"
          "      <prefix>.scmask       pmemd scmask1 / scmask2\n"
          "      <prefix>.atoms        per-slot kind, names, charges\n"
          "    series also: <outprefix><name>.sorted.{lib|mol2|pdb},\n"
          "                 <tioutprefix><name>.*, <mapoutprefix><name>.map,\n"
          "                 parentout score table.\n"
          "\n"
          "  Matching: seed auto|names|na|aa|none (default auto). Insertions bonded\n"
          "  to a mapped parent follow that parent; other leftovers go before\n"
          "  anchor <atomname> (default O3').\n"
          "\n"
          "  Inputs: Amber OFF .lib (readdata → Name[Unit]), or mol2 / pdb / Amber\n"
          "  topology via parm — or pass a file path and timap loads it. Formats may\n"
          "  be mixed; everything is matched as Topology. ModXNA parents and\n"
          "  ff19SB amino19.lib: $CPPTRAJHOME/dat/timap/. Immediate command; add go\n"
          "  if the input has trajin. Alias: timatch.\n"
          "\n"
          "  Examples:\n"
          "    > readdata parent.lib name parent\n"
          "    > readdata analog.lib name analog\n"
          "    > timap analog[analog] template parent[parent] out analog.lib\n"
          "    > timap analog[analog] naorder out analog.lib\n"
          "    > timap ncaa[ncaa] aaorder out ncaa.lib\n"
          "    > timap analog[analog] template parent[parent] tiout analog_ti\n"
          "    > parm a.mol2 name A\n"
          "    > parm b.pdb name B\n"
          "    > timap series A B outprefix s_ tioutprefix s_ti_ parentout parent.dat\n"
          "    > timap phenol.mol2 template benzene.pdb out phenol.lib\n"
          "    > timap phenol.mol2 template benzene.pdb out phenol.mol2\n"
          "    > timap phenol.mol2 template benzene.pdb out phenol.pdb\n"
          "    > timap phenol.mol2 template cys.lib[CYS] mapout mix.map maponly\n");
}

/** Look up a topology by the name the user typed (already-loaded sets only).
  *
  * Order of attempts:
  *   1. A TOPOLOGY dataset (parm / parmwrite sets: mol2, pdb, Amber, ...).
  *   2. A COORDS dataset — Amber OFF units loaded with readdata appear as
  *      LibName[UnitName], e.g. FLE[FLE] or CYS[CYS].
  *   3. A purely numeric string, treated as a parm index.
  *
  * \param dsOut If non-null, receives the dataset pointer (TOPOLOGY or COORDS).
  */
static Topology* FindNamedTop(DataSetList& dsl, std::string const& name, DataSet** dsOut)
{
  if (dsOut != 0) *dsOut = 0;
  if (name.empty()) return 0;
  DataSet* ds = dsl.FindSetOfType(name, DataSet::TOPOLOGY);
  if (ds != 0) {
    if (dsOut != 0) *dsOut = ds;
    return ((DataSet_Topology*)ds)->TopPtr();
  }
  ds = dsl.FindSetOfGroup(name, DataSet::COORDINATES);
  if (ds != 0) {
    if (dsOut != 0) *dsOut = ds;
    return ((DataSet_Coords*)ds)->TopPtr();
  }
  bool isInt = !name.empty();
  for (size_t i = 0; i < name.size(); i++) {
    if (name[i] < '0' || name[i] > '9') { isInt = false; break; }
  }
  if (!isInt) return 0;
  ArgList tmp(name);
  return dsl.GetTopByIndex(tmp);
}

/** Split optional Amber-OFF unit bracket: path[Unit] → path, Unit. */
static void ParseTopSpec(std::string const& spec, std::string& primary, std::string& unit)
{
  primary = spec;
  unit.clear();
  if (spec.empty()) return;
  size_t rb = spec.rfind(']');
  size_t lb = spec.rfind('[');
  if (lb != std::string::npos && rb != std::string::npos &&
      rb == spec.size() - 1 && rb > lb)
  {
    primary = spec.substr(0, lb);
    unit = spec.substr(lb + 1, rb - lb - 1);
  }
}

/** Basename without extension (phenol.mol2 → phenol). */
static std::string FileStem(FileName const& fn)
{
  std::string b = fn.Base();
  std::string e = fn.Ext();
  if (!e.empty() && b.size() > e.size() &&
      b.compare(b.size() - e.size(), e.size(), e) == 0)
    return b.substr(0, b.size() - e.size());
  return b;
}

/**
 * Resolve a user topology specifier to Topology.
 * Accepts an already-loaded set name, or a mol2/pdb/Amber/.lib file path
 * (auto-loaded via ParmFile or DataFile Amber OFF). All formats become the
 * same internal Topology for TIMatch; mix-and-match is allowed.
 */
static Topology* ResolveTop(CpptrajState& State, std::string const& spec, DataSet** dsOut)
{
  if (dsOut != 0) *dsOut = 0;
  if (spec.empty()) return 0;

  Topology* top = FindNamedTop(State.DSL(), spec, dsOut);
  if (top != 0) return top;

  std::string primary, unit;
  ParseTopSpec(spec, primary, unit);
  if (primary != spec) {
    top = FindNamedTop(State.DSL(), primary, dsOut);
    if (top != 0) return top;
    if (!unit.empty()) {
      std::string withUnit = primary + "[" + unit + "]";
      top = FindNamedTop(State.DSL(), withUnit, dsOut);
      if (top != 0) return top;
    }
  }

  bool maybeFile = File::Exists(primary) ||
                   primary.find('.') != std::string::npos ||
                   primary.find('/') != std::string::npos ||
                   primary.find('\\') != std::string::npos;
  if (!maybeFile || !File::Exists(primary))
    return 0;

  FileName fn(primary);
  std::string stem = FileStem(fn);
  if (stem.empty()) stem = fn.Base();

  // Topology formats: mol2, pdb, Amber parm, PSF, ...
  if (ParmFile::DetectFormat(fn) != ParmFile::UNKNOWN_PARM) {
    ArgList nameArg(std::string("name ") + stem);
    if (State.AddTopology(primary, nameArg)) {
      mprinterr("Error: timap: could not load topology file '%s'\n", primary.c_str());
      return 0;
    }
    mprintf("\tLoaded topology file '%s' as '%s'\n", primary.c_str(), stem.c_str());
    return FindNamedTop(State.DSL(), stem, dsOut);
  }

  // Amber OFF .lib (and other DataFile formats that yield COORDS with Topology)
  DataFile dfile;
  ArgList nameArg(std::string("name ") + stem);
  if (dfile.ReadDataIn(fn, nameArg, State.DSL())) {
    mprinterr("Error: timap: could not load data/topology file '%s'\n", primary.c_str());
    return 0;
  }
  std::string lookup = unit.empty() ? stem : (stem + "[" + unit + "]");
  top = FindNamedTop(State.DSL(), lookup, dsOut);
  if (top == 0 && !unit.empty())
    top = FindNamedTop(State.DSL(), stem, dsOut);
  if (top == 0) {
    mprinterr("Error: timap: loaded '%s' but found no topology set named '%s'\n",
              primary.c_str(), lookup.c_str());
    return 0;
  }
  mprintf("\tLoaded file '%s' as '%s'\n", primary.c_str(), lookup.c_str());
  return top;
}

/** Fill frm with coordinates for top. Prefer a COORDS frame; else the original file. */
static bool LoadCoords(Topology const& top, DataSet* ds, Frame& frm)
{
  frm.SetupFrame(top.Natom());
  for (int i = 0; i < frm.Natom() * 3; i++)
    frm[i] = 0.0;
  if (ds != 0 && ds->Group() == DataSet::COORDINATES && ds->Size() > 0) {
    DataSet_Coords* crd = (DataSet_Coords*)ds;
    Frame tmp = crd->AllocateFrame();
    crd->GetFrame(0, tmp);
    int ncopy = std::min(tmp.Natom(), frm.Natom());
    for (int i = 0; i < ncopy * 3; i++)
      frm[i] = tmp[i];
    return ncopy == frm.Natom();
  }
  if (top.OriginalFilename().empty())
    return frm.Natom() == 0;
  Trajin_Single traj;
  ArgList empty;
  Topology* tptr = const_cast<Topology*>(&top);
  if (traj.SetupTrajRead(top.OriginalFilename(), empty, tptr))
    return false;
  if (traj.BeginTraj())
    return false;
  int got = traj.GetNextFrame(frm);
  traj.EndTraj();
  return got == 1;
}

static void CopyXyz(Frame const& src, int oldat, Frame& dst, int newat)
{
  if (oldat < 0 || newat < 0) return;
  if (oldat >= src.Natom() || newat >= dst.Natom()) return;
  dst[newat * 3    ] = src[oldat * 3    ];
  dst[newat * 3 + 1] = src[oldat * 3 + 1];
  dst[newat * 3 + 2] = src[oldat * 3 + 2];
}

/** Permute src coordinates into dst using Map[new] = old. */
static void ApplyOrderToFrame(Frame const& src, std::vector<int> const& order, Frame& dst)
{
  dst.SetupFrame((int)order.size());
  for (int i = 0; i < (int)order.size() * 3; i++)
    dst[i] = 0.0;
  for (int i = 0; i < (int)order.size(); i++)
    CopyXyz(src, order[i], dst, i);
}

/** Amber OFF unit names: letters, digits, underscore only. */
static std::string OffUnitName(std::string in, char fallback)
{
  std::string out;
  for (size_t i = 0; i < in.size(); i++) {
    char c = in[i];
    if (std::isalnum((unsigned char)c) || c == '_')
      out += c;
  }
  if (out.empty()) {
    out = "L";
    out += fallback;
  }
  if (out.size() > 20) out.resize(20);
  return out;
}

static NameType ResNameOf(Topology const& top, const char* fallback)
{
  if (top.Nres() > 0)
    return top.Res(0).Name();
  return NameType(fallback);
}

/** Build one TI end-state topology+frame from the dual-slot list.
  * lambda0: template atoms real, target insertions dummy (charge 0, mass 0, type DUM).
  * lambda1: target atoms real, unmatched template atoms dummy.
  */
static int BuildTiUnit(bool lambda0,
                       Topology const& tgt, Topology const& tpl,
                       Frame const& tgtX, Frame const& tplX,
                       std::vector<TIMatch::Result::DualSlot> const& dual,
                       Topology& outTop, Frame& outFrm,
                       NameType const& resname)
{
  typedef TIMatch::Result::DualSlot Slot;
  outTop = Topology();
  Residue res(resname, 1, ' ', "");
  int nd = (int)dual.size();
  outFrm.SetupFrame(nd);
  for (int i = 0; i < nd * 3; i++)
    outFrm[i] = 0.0;
  std::vector<int> tgtDual(tgt.Natom(), -1);
  std::vector<int> tplDual(tpl.Natom(), -1);
  for (int i = 0; i < nd; i++) {
    Slot const& s = dual[i];
    if (s.tgt_ >= 0) tgtDual[s.tgt_] = i;
    if (s.tpl_ >= 0) tplDual[s.tpl_] = i;
    Atom src;
    bool dummy = false;
    int xyzFrom = -1;
    Frame const* xyzSrc = 0;
    if (lambda0) {
      if (s.tpl_ >= 0) {
        src = tpl[s.tpl_];
        xyzSrc = &tplX;
        xyzFrom = s.tpl_;
      } else {
        src = tgt[s.tgt_];
        xyzSrc = &tgtX;
        xyzFrom = s.tgt_;
        dummy = true;
      }
    } else {
      if (s.tgt_ >= 0) {
        src = tgt[s.tgt_];
        xyzSrc = &tgtX;
        xyzFrom = s.tgt_;
      } else {
        src = tpl[s.tpl_];
        xyzSrc = &tplX;
        xyzFrom = s.tpl_;
        dummy = true;
      }
    }
    src.ClearBonds();
    if (dummy) {
      src.SetCharge(0.0);
      src.SetMass(0.0);
      src.SetTypeName("DUM");
    }
    outTop.AddTopAtom(src, res);
    if (xyzSrc != 0)
      CopyXyz(*xyzSrc, xyzFrom, outFrm, i);
  }

  if (lambda0) {
    for (int r = 0; r < tpl.Natom(); r++) {
      if (tplDual[r] < 0) continue;
      Atom const& at = tpl[r];
      for (Atom::bond_iterator b = at.bondbegin(); b != at.bondend(); ++b) {
        if (*b > r && tplDual[*b] >= 0)
          outTop.AddBond(tplDual[r], tplDual[*b], -1);
      }
    }
    for (int i = 0; i < nd; i++) {
      if (!dual[i].IsTgtOnly()) continue;
      int n = dual[i].tgt_;
      int attach = -1;
      Atom const& at = tgt[n];
      for (Atom::bond_iterator b = at.bondbegin(); b != at.bondend(); ++b) {
        int d = tgtDual[*b];
        if (d < 0) continue;
        if (attach < 0 || dual[d].IsShared())
          attach = d;
      }
      if (attach >= 0)
        outTop.AddBond(i, attach, -1);
      else
        mprintf("Warning: timap: dummy %s has no bonded neighbor in the dual layout.\n",
                outTop[i].c_str());
    }
  } else {
    for (int n = 0; n < tgt.Natom(); n++) {
      if (tgtDual[n] < 0) continue;
      Atom const& at = tgt[n];
      for (Atom::bond_iterator b = at.bondbegin(); b != at.bondend(); ++b) {
        if (*b > n && tgtDual[*b] >= 0)
          outTop.AddBond(tgtDual[n], tgtDual[*b], -1);
      }
    }
    for (int i = 0; i < nd; i++) {
      if (!dual[i].IsTplOnly()) continue;
      int r = dual[i].tpl_;
      int attach = -1;
      Atom const& at = tpl[r];
      for (Atom::bond_iterator b = at.bondbegin(); b != at.bondend(); ++b) {
        int d = tplDual[*b];
        if (d < 0) continue;
        if (attach < 0 || dual[d].IsShared())
          attach = d;
      }
      if (attach >= 0)
        outTop.AddBond(i, attach, -1);
      else
        mprintf("Warning: timap: dummy %s has no bonded neighbor in the dual layout.\n",
                outTop[i].c_str());
    }
  }
  outTop.CommonSetup();
  return 0;
}

static int WriteMol2File(std::string const& fname, Topology& top, Frame const& frm,
                         DataSetList const& dsl)
{
  Trajout_Single out;
  ArgList empty;
  if (out.PrepareTrajWrite(fname, empty, dsl, &top, CoordinateInfo(), 1,
                           TrajectoryFile::MOL2FILE))
  {
    mprinterr("Error: timap: could not set up mol2 '%s'\n", fname.c_str());
    return 1;
  }
  if (out.WriteSingle(0, frm)) {
    mprinterr("Error: timap: writing mol2 '%s'\n", fname.c_str());
    return 1;
  }
  out.EndTraj();
  mprintf("\tWrote '%s' (%i atoms).\n", fname.c_str(), top.Natom());
  return 0;
}

/** Write PDB via Traj_PDBfile; include CONECT so bonds round-trip for rematching. */
static int WritePdbFile(std::string const& fname, Topology& top, Frame const& frm,
                        DataSetList const& dsl)
{
  Trajout_Single out;
  ArgList pdbArgs("conect");
  if (out.PrepareTrajWrite(fname, pdbArgs, dsl, &top, CoordinateInfo(), 1,
                           TrajectoryFile::PDBFILE))
  {
    mprinterr("Error: timap: could not set up pdb '%s'\n", fname.c_str());
    return 1;
  }
  if (out.WriteSingle(0, frm)) {
    mprinterr("Error: timap: writing pdb '%s'\n", fname.c_str());
    return 1;
  }
  out.EndTraj();
  mprintf("\tWrote '%s' (%i atoms).\n", fname.c_str(), top.Natom());
  return 0;
}

/// Aligned-structure output format (default Amber OFF .lib for LEaP/TI).
enum AlignedOutFmt { AOUT_LIB = 0, AOUT_MOL2, AOUT_PDB };

static const char* AlignedOutExt(AlignedOutFmt f)
{
  switch (f) {
    case AOUT_MOL2: return ".sorted.mol2";
    case AOUT_PDB:  return ".sorted.pdb";
    case AOUT_LIB:
    default:        return ".sorted.lib";
  }
}

static const char* AlignedOutLabel(AlignedOutFmt f)
{
  switch (f) {
    case AOUT_MOL2: return "mol2";
    case AOUT_PDB:  return "pdb";
    case AOUT_LIB:
    default:        return "Amber OFF .lib";
  }
}

static int ParseAlignedOutFmt(std::string const& s, AlignedOutFmt& out)
{
  std::string l = ToLower(s);
  if (l == "lib" || l == "off" || l == "amberlib") { out = AOUT_LIB;  return 0; }
  if (l == "mol2")                                   { out = AOUT_MOL2; return 0; }
  if (l == "pdb")                                    { out = AOUT_PDB;  return 0; }
  return 1;
}

/** Prefer extension of out <file>; otherwise keep def (from fmt keyword). */
static AlignedOutFmt FmtFromFilename(std::string const& fname, AlignedOutFmt def)
{
  FileName fn(fname);
  std::string e = ToLower(fn.Ext());
  if (e == ".mol2") return AOUT_MOL2;
  if (e == ".pdb")  return AOUT_PDB;
  if (e == ".lib" || e == ".off") return AOUT_LIB;
  return def;
}

static std::string Q(std::string const& s)
{
  return std::string("\"") + s + "\"";
}

static std::string AtomTypeStr(Atom const& a)
{
  std::string t = a.Type().Truncated();
  if (t.empty() && a.ElementName() != 0)
    t = a.ElementName();
  return t;
}

/** 1-based index of the first atom with this name, or 0 if missing. */
static int AtomNum1(Topology const& top, const char* n)
{
  for (int i = 0; i < top.Natom(); i++) {
    if (top[i].Name() == n)
      return i + 1;
  }
  return 0;
}

/** Head/tail connect atoms for LEaP (P/O3' or N/C). */
static void LeapConnect(Topology const& top, int& c1, int& c2, const char*& restype)
{
  c1 = AtomNum1(top, "P");
  c2 = AtomNum1(top, "O3'");
  if (c1 > 0 && c2 > 0) {
    restype = "n";
    return;
  }
  c1 = AtomNum1(top, "N");
  c2 = AtomNum1(top, "C");
  if (c1 > 0 && c2 > 0) {
    restype = "p";
    return;
  }
  c1 = 0;
  c2 = 0;
  restype = "?";
}

/** Minimal Amber OFF writer for a single-residue TI unit. */
static int WriteAmberLib(std::string const& fname, std::string const& unit,
                         Topology const& top, Frame const& frm)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open '%s' for Amber OFF write.\n", fname.c_str());
    return 1;
  }
  int n = top.Natom();
  std::string res = (top.Nres() > 0) ? top.Res(0).Name().Truncated() : unit;
  int c1, c2;
  const char* restype;
  LeapConnect(top, c1, c2, restype);
  out.Printf("!!index array str\n \"%s\"\n", unit.c_str());
  out.Printf("!entry.%s.unit.atoms table  str name  str type  int typex  "
             "int resx  int flags  int seq  int elmnt  dbl chg\n", unit.c_str());
  for (int i = 0; i < n; i++) {
    Atom const& a = top[i];
    out.Printf(" %s %s %5d %3d %7d %3d %3d %10.6f\n",
               Q(a.Name().Truncated()).c_str(),
               Q(AtomTypeStr(a)).c_str(),
               0, 1, 131072, i + 1, a.AtomicNumber(), a.Charge());
  }
  out.Printf("!entry.%s.unit.atomspertinfo table  str pname  str ptype  "
             "int ptypex  int pelmnt  dbl pchg\n", unit.c_str());
  for (int i = 0; i < n; i++)
    out.Printf(" %s %s 0 -1 0.0\n",
               Q(top[i].Name().Truncated()).c_str(),
               Q(AtomTypeStr(top[i])).c_str());
  out.Printf("!entry.%s.unit.boundbox array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n",
             unit.c_str());
  out.Printf("!entry.%s.unit.childsequence single int\n 2\n", unit.c_str());
  out.Printf("!entry.%s.unit.connect array int\n %d\n %d\n", unit.c_str(), c1, c2);
  out.Printf("!entry.%s.unit.connectivity table  int atom1x  int atom2x  int flags\n",
             unit.c_str());
  for (int i = 0; i < n; i++) {
    Atom const& a = top[i];
    for (Atom::bond_iterator b = a.bondbegin(); b != a.bondend(); ++b) {
      if (*b > i)
        out.Printf(" %d %d 1\n", i + 1, (*b) + 1);
    }
  }
  out.Printf("!entry.%s.unit.hierarchy table  str abovetype  int abovex  "
             "str belowtype  int belowx\n", unit.c_str());
  out.Printf(" \"U\" 0 \"R\" 1\n");
  for (int i = 0; i < n; i++)
    out.Printf(" \"R\" 1 \"A\" %d\n", i + 1);
  out.Printf("!entry.%s.unit.name single str\n %s\n", unit.c_str(), Q(unit).c_str());
  out.Printf("!entry.%s.unit.positions table  dbl x  dbl y  dbl z\n", unit.c_str());
  for (int i = 0; i < n; i++)
    out.Printf(" %.6f %.6f %.6f\n", frm[i * 3], frm[i * 3 + 1], frm[i * 3 + 2]);
  out.Printf("!entry.%s.unit.residueconnect table  int c1x  int c2x  int c3x  "
             "int c4x  int c5x  int c6x\n %d %d 0 0 0 0\n", unit.c_str(), c1, c2);
  out.Printf("!entry.%s.unit.residues table  str name  int seq  int childseq  "
             "int startatomx  str restype  int imagingx\n", unit.c_str());
  out.Printf(" %s 1 %d 1 \"%s\" 0\n", Q(res).c_str(), n + 1, restype);
  out.Printf("!entry.%s.unit.residuesPdbSequenceNumber array int\n 0\n", unit.c_str());
  out.Printf("!entry.%s.unit.solventcap array dbl\n -1.000000\n 0.0\n 0.0\n 0.0\n 0.0\n",
             unit.c_str());
  out.Printf("!entry.%s.unit.velocities table  dbl x  dbl y  dbl z\n", unit.c_str());
  for (int i = 0; i < n; i++)
    out.Printf(" 0.0 0.0 0.0\n");
  out.CloseFile();
  mprintf("\tWrote '%s' (unit %s, %i atoms).\n", fname.c_str(), unit.c_str(), n);
  return 0;
}

/** Write the remapped structure in the requested format (lib / mol2 / pdb). */
static int WriteAlignedOut(std::string const& fname, AlignedOutFmt fmt,
                           Topology& top, Frame const& frm, DataSetList const& dsl)
{
  std::string unit = OffUnitName(ResNameOf(top, "TGT").Truncated(), 'T');
  top.SetParmName(unit, FileName(fname));
  switch (fmt) {
    case AOUT_MOL2:
      return WriteMol2File(fname, top, frm, dsl);
    case AOUT_PDB:
      return WritePdbFile(fname, top, frm, dsl);
    case AOUT_LIB:
    default:
      return WriteAmberLib(fname, unit, top, frm);
  }
}

static int WriteScmask(std::string const& fname,
                       Topology const& top0, Topology const& top1,
                       std::vector<TIMatch::Result::DualSlot> const& dual)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  int nShared = 0, nD0 = 0, nD1 = 0;
  for (size_t i = 0; i < dual.size(); i++) {
    if (dual[i].IsShared()) nShared++;
    else if (dual[i].IsTgtOnly()) nD0++;
    else nD1++;
  }
  out.Printf("# timap dual-topology TI masks\n");
  out.Printf("# n_dual= %zu  n_shared= %i  dummy_in_lambda0= %i  dummy_in_lambda1= %i\n",
             dual.size(), nShared, nD0, nD1);
  out.Printf("# lambda 0 = template (real unmatched, dummy insertions)\n");
  out.Printf("# lambda 1 = target   (real insertions, dummy unmatched)\n");
  out.Printf("# pmemd: scmask1 is unique to lambda 0; scmask2 is unique to lambda 1.\n");
  std::string m1, m2;
  for (int i = 0; i < (int)dual.size(); i++) {
    if (dual[i].IsTplOnly()) {
      if (!m1.empty()) m1 += ",";
      m1 += ":1@";
      m1 += top0[i].Name().Truncated();
    }
    if (dual[i].IsTgtOnly()) {
      if (!m2.empty()) m2 += ",";
      m2 += ":1@";
      m2 += top1[i].Name().Truncated();
    }
  }
  out.Printf("scmask1=\"%s\"\n", m1.c_str());
  out.Printf("scmask2=\"%s\"\n", m2.c_str());
  out.CloseFile();
  mprintf("\tWrote '%s'\n", fname.c_str());
  mprintf("\tscmask1 (lambda 0 unique) = %s\n", m1.empty() ? "(none)" : m1.c_str());
  mprintf("\tscmask2 (lambda 1 unique) = %s\n", m2.empty() ? "(none)" : m2.c_str());
  return 0;
}

/** Per-slot dual-topology table: kind, names, and charges in both end states. */
static int WriteDualAtoms(std::string const& fname,
                          Topology const& top0, Topology const& top1,
                          std::vector<TIMatch::Result::DualSlot> const& dual)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  out.Printf("# timap dual-topology atoms\n");
  out.Printf("# Kind: SHARED = real in both; TPL_ONLY = dummy in lambda 1; "
             "TGT_ONLY = dummy in lambda 0\n");
  out.Printf("%-4s %-8s %-8s %10s %-8s %10s\n",
             "#Idx", "Kind", "Name0", "Q0", "Name1", "Q1");
  for (int i = 0; i < (int)dual.size(); i++) {
    const char* kind = "SHARED";
    if (dual[i].IsTplOnly()) kind = "TPL_ONLY";
    else if (dual[i].IsTgtOnly()) kind = "TGT_ONLY";
    out.Printf(" %3i %-8s %-8s %10.6f %-8s %10.6f\n",
               i + 1, kind,
               top0[i].Name().Truncated().c_str(), top0[i].Charge(),
               top1[i].Name().Truncated().c_str(), top1[i].Charge());
  }
  out.CloseFile();
  mprintf("\tWrote '%s'\n", fname.c_str());
  return 0;
}

/** Write the human-readable correspondence used by Test_TIMap. */
static int WriteMapFile(std::string const& fname, Topology const& tgt, Topology const& tpl,
                        TIMatch::Result const& R, std::string const& tgtName,
                        std::string const& tplName)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open map file '%s'\n", fname.c_str());
    return 1;
  }
  out.Printf("# timap tgt='%s' template='%s'\n", tgtName.c_str(), tplName.c_str());
  out.Printf("# kind tgt=%s template=%s\n", R.tgtKind_.c_str(), R.tplKind_.c_str());
  out.Printf("# mapped= %i  insertion= %i  unmapped_template= %i  n_tgt= %i  n_tpl= %i\n",
             R.nMapped_, R.nInsertion_, R.nUnmappedTpl_, tgt.Natom(), tpl.Natom());
  int aWidth = std::max(6, DigitWidth(tgt.Natom()));
  aWidth = std::max(aWidth, DigitWidth(tpl.Natom()));
  int nWidth = 6;
  for (int i = 0; i < tgt.Natom(); i++)
    nWidth = std::max(nWidth, tgt[i].Name().len());
  for (int i = 0; i < tpl.Natom(); i++)
    nWidth = std::max(nWidth, tpl[i].Name().len());
  out.Printf("%-*s %*s %-*s %*s %-*s %s\n",
             aWidth, "#Out",
             aWidth, "TgtAt",
             nWidth, "TgtName",
             aWidth, "TplAt",
             nWidth, "TplName",
             "Ins");
  for (int k = 0; k < (int)R.outputOrder_.size(); k++) {
    int oldat = R.outputOrder_[k];
    int tplidx = R.mapping_[oldat];
    if (tplidx < 0)
      out.Printf("%*i %*i %-*s %*s %-*s %s\n",
                 aWidth, k + 1,
                 aWidth, oldat + 1,
                 nWidth, tgt[oldat].c_str(),
                 aWidth, "---",
                 nWidth, "---",
                 "1");
    else
      out.Printf("%*i %*i %-*s %*i %-*s %s\n",
                 aWidth, k + 1,
                 aWidth, oldat + 1,
                 nWidth, tgt[oldat].c_str(),
                 aWidth, tplidx + 1,
                 nWidth, tpl[tplidx].c_str(),
                 "0");
  }
  if (R.nUnmappedTpl_ > 0) {
    out.Printf("# unmapped template atoms:\n");
    std::vector<char> used(tpl.Natom(), 0);
    for (int n = 0; n < tgt.Natom(); n++) {
      if (R.mapping_[n] >= 0) used[R.mapping_[n]] = 1;
    }
    for (int r = 0; r < tpl.Natom(); r++) {
      if (!used[r])
        out.Printf("#   %i %s\n", r + 1, tpl[r].c_str());
    }
  }
  out.CloseFile();
  return 0;
}

/** Filename-safe token from a dataset name (FLE[FLE] -> FLE_FLE). */
static std::string SafeFileToken(std::string const& in)
{
  std::string out;
  for (size_t i = 0; i < in.size(); i++) {
    char c = in[i];
    if (std::isalnum((unsigned char)c) || c == '_' || c == '-' || c == '.')
      out += c;
    else
      out += '_';
  }
  if (out.empty()) out = "mol";
  return out;
}

/** Append comma-split tokens from one user string into names. */
static void PushSeriesName(std::vector<std::string>& names, std::string const& raw)
{
  if (raw.empty()) return;
  if (raw.find(',') == std::string::npos) {
    names.push_back(raw);
    return;
  }
  ArgList parts(raw, ",");
  for (int i = 0; i < parts.Nargs(); i++) {
    std::string s = parts[i];
    // Trim spaces
    size_t b = 0;
    while (b < s.size() && (s[b] == ' ' || s[b] == '\t')) b++;
    size_t e = s.size();
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t')) e--;
    if (e > b) names.push_back(s.substr(b, e - b));
  }
}

/**
 * leadopt-style pairwise similarity in [0,1].
 * Adapted from util/leadopt/similarity.py::exp_delta / by_heavy_atom_count
 * (Summa Lab / Redmann & Summa): score = exp(-BETA * delta), BETA = 0.1,
 * where delta counts atoms not in the common substructure. Here delta is
 * TIMatch insertions + unmatched template atoms (partial map proxy for MCS).
 */
static double LeadoptSimScore(int nInsertion, int nUnmappedTpl)
{
  const double BETA = 0.1;
  return std::exp(-BETA * (double)(nInsertion + nUnmappedTpl));
}

struct SeriesMember {
  std::string name;
  Topology* top;
  DataSet* ds;
};

/**
 * Choose the series parent: maximize total mapped atoms when every other
 * member is mapped onto the candidate ("shares the most"); on ties prefer
 * the minimum structure (fewest atoms), then the leadopt-style similarity sum.
 * Inspired by util/leadopt pairwise MCS scoring (mcs.py / similarity.py /
 * graph.py) and the star-parent use case of GraphGenerator4.py (Redmann &
 * Summa, U. New Orleans) without importing that TI-path planner.
 * \return parent index, or -1 on failure.
 */
static int SelectSeriesParent(TIMatch const& matcher,
                              std::vector<SeriesMember> const& mem,
                              std::vector<double>& outSimSum,
                              std::vector<int>& outMappedSum)
{
  int n = (int)mem.size();
  outSimSum.assign(n, 0.0);
  outMappedSum.assign(n, 0);
  if (n < 2) return -1;
  int best = -1;
  for (int i = 0; i < n; i++) {
    double simSum = 0.0;
    int mappedSum = 0;
    for (int j = 0; j < n; j++) {
      if (j == i) continue;
      TIMatch::Result R;
      if (matcher.Match(*mem[j].top, *mem[i].top, R)) {
        mprinterr("Error: timap: series match of '%s' onto candidate '%s' failed.\n",
                  mem[j].name.c_str(), mem[i].name.c_str());
        return -1;
      }
      mappedSum += R.nMapped_;
      simSum += LeadoptSimScore(R.nInsertion_, R.nUnmappedTpl_);
    }
    outSimSum[i] = simSum;
    outMappedSum[i] = mappedSum;
    if (best < 0) {
      best = i;
      continue;
    }
    bool better = false;
    if (mappedSum > outMappedSum[best])
      better = true;
    else if (mappedSum == outMappedSum[best]) {
      if (mem[i].top->Natom() < mem[best].top->Natom())
        better = true;
      else if (mem[i].top->Natom() == mem[best].top->Natom() &&
               simSum > outSimSum[best] + 1.0e-12)
        better = true;
    }
    if (better) best = i;
  }
  return best;
}

/** Parse arguments, run TIMatch::Match, write an aligned structure (lib/mol2/pdb),
  * and optionally dual-topology TI files. series mode auto-selects a parent.
  */
// Exec_TIMap::Execute()
Exec::RetType Exec_TIMap::Execute(CpptrajState& State, ArgList& argIn) {
  std::string mapout = argIn.GetStringKey("mapout");
  std::string libout = argIn.GetStringKey("out");
  std::string newname = argIn.GetStringKey("name");
  std::string tplName = argIn.GetStringKey("template");
  std::string seedStr = argIn.GetStringKey("seed");
  std::string anchor = argIn.GetStringKey("anchor");
  std::string tiout = argIn.GetStringKey("tiout");
  std::string outprefix = argIn.GetStringKey("outprefix");
  std::string tioutprefix = argIn.GetStringKey("tioutprefix");
  std::string mapoutprefix = argIn.GetStringKey("mapoutprefix");
  std::string parentout = argIn.GetStringKey("parentout");
  std::string fmtStr = argIn.GetStringKey("fmt");
  bool maponly = argIn.hasKey("maponly");
  bool replace = argIn.hasKey("replace");
  bool naorder = argIn.hasKey("naorder");
  bool aaorder = argIn.hasKey("aaorder");
  bool series = argIn.hasKey("series");
  if (naorder && aaorder) {
    mprinterr("Error: timap: specify naorder or aaorder, not both.\n");
    return CpptrajState::ERR;
  }
  AlignedOutFmt outFmt = AOUT_LIB;
  if (!fmtStr.empty() && ParseAlignedOutFmt(fmtStr, outFmt)) {
    mprinterr("Error: timap: unrecognized fmt '%s' (lib|mol2|pdb).\n", fmtStr.c_str());
    return CpptrajState::ERR;
  }
  if (!libout.empty())
    outFmt = FmtFromFilename(libout, outFmt);

  // Collect remaining unmarked names (series members, or tgt [/ template]).
  std::vector<std::string> names;
  std::string next = argIn.GetStringNext();
  while (!next.empty()) {
    PushSeriesName(names, next);
    next = argIn.GetStringNext();
  }
  // Legacy: template <name> already consumed; if no series and template empty,
  // second positional was historically the optional template without the keyword.
  // That form is: timap <tgt> <tpl> — already collected into names.

  if (series) {
    if (!tplName.empty()) {
      mprinterr("Error: timap: series chooses the parent; do not also give template.\n");
      return CpptrajState::ERR;
    }
    if (replace || !newname.empty()) {
      mprinterr("Error: timap: series does not support replace / name.\n");
      return CpptrajState::ERR;
    }
    if (names.size() < 2) {
      mprinterr("Error: timap: series needs at least two topologies.\n");
      return CpptrajState::ERR;
    }

    TIMatch matcher;
    matcher.SetDebug(State.Debug());
    matcher.SetUseNaOrder(naorder);
    matcher.SetUseAaOrder(aaorder);
    if (aaorder && anchor.empty())
      matcher.SetAnchorName("C");
    if (!anchor.empty()) matcher.SetAnchorName(anchor);
    if (!seedStr.empty()) {
      std::string l = ToLower(seedStr);
      if (l != "auto" && l != "names" && l != "na" && l != "aa" && l != "none") {
        mprinterr("Error: timap: unrecognized seed '%s'\n", seedStr.c_str());
        return CpptrajState::ERR;
      }
      matcher.SetSeed(TIMatch::SeedFromString(seedStr));
    }

    std::vector<SeriesMember> mem;
    mem.reserve(names.size());
    for (size_t i = 0; i < names.size(); i++) {
      SeriesMember m;
      m.name = names[i];
      m.ds = 0;
      m.top = ResolveTop(State, m.name, &m.ds);
      if (m.top == 0) {
        mprinterr("Error: timap: series member '%s' not found.\n", m.name.c_str());
        return CpptrajState::ERR;
      }
      mem.push_back(m);
    }

    mprintf("    TIMAP: Series of %zu topologies; selecting parent.\n", mem.size());
    mprintf("\tSeed: %s\n", TIMatch::SeedStr(
              seedStr.empty() ? TIMatch::SEED_AUTO
                              : TIMatch::SeedFromString(seedStr)));
    if (naorder)
      mprintf("\tOrder: nucleic-acid canonical walk of the selected parent.\n");
    else if (aaorder)
      mprintf("\tOrder: amino-acid canonical walk (ff19SB) of the selected parent.\n");
    else
      mprintf("\tOrder: selected parent's current atom order.\n");
    mprintf("\tParent score: maximize total mapped atoms; ties prefer fewer atoms,\n"
            "\t  then leadopt-style exp(-0.1*(ins+unmap)) (see util/leadopt citation).\n");

    std::vector<double> simSum;
    std::vector<int> mappedSum;
    int pIdx = SelectSeriesParent(matcher, mem, simSum, mappedSum);
    if (pIdx < 0) return CpptrajState::ERR;

    SeriesMember const& parent = mem[pIdx];
    mprintf("\tSelected parent: '%s' (%i atoms, sim-sum=%.6f, mapped-sum=%i).\n",
            parent.name.c_str(), parent.top->Natom(), simSum[pIdx], mappedSum[pIdx]);
    for (size_t i = 0; i < mem.size(); i++) {
      mprintf("\t  candidate %-16s atoms=%4i  sim-sum=%.6f  mapped-sum=%i%s\n",
              mem[i].name.c_str(), mem[i].top->Natom(),
              simSum[i], mappedSum[i],
              ((int)i == pIdx) ? "  <-- parent" : "");
    }

    if (!parentout.empty()) {
      CpptrajFile pout;
      if (pout.OpenWrite(parentout)) {
        mprinterr("Error: timap: could not write parentout '%s'\n", parentout.c_str());
        return CpptrajState::ERR;
      }
      pout.Printf("# timap series parent selection\n");
      pout.Printf("# score = sum_j exp(-0.1*(insertions+unmapped)) mapping j onto candidate\n");
      pout.Printf("# adapted from util/leadopt/similarity.py (Redmann & Summa / Summa Lab)\n");
      pout.Printf("parent %s\n", parent.name.c_str());
      pout.Printf("#Idx Name                 Natom   SimSum  MappedSum\n");
      for (size_t i = 0; i < mem.size(); i++) {
        pout.Printf(" %3zu %-20s %5i %10.6f %10i%s\n",
                    i + 1, mem[i].name.c_str(), mem[i].top->Natom(),
                    simSum[i], mappedSum[i],
                    ((int)i == pIdx) ? " parent" : "");
      }
      pout.CloseFile();
      mprintf("\tWrote parent selection to '%s'\n", parentout.c_str());
    }

    bool writeLib = !maponly;
    for (size_t i = 0; i < mem.size(); i++) {
      SeriesMember const& tgt = mem[i];
      TIMatch::Result R;
      if (matcher.Match(*tgt.top, *parent.top, R)) return CpptrajState::ERR;

      mprintf("\t[%zu/%zu] '%s' onto parent '%s': mapped %i / %i "
              "(%i insertions, %i parent unmatched).\n",
              i + 1, mem.size(), tgt.name.c_str(), parent.name.c_str(),
              R.nMapped_, tgt.top->Natom(), R.nInsertion_, R.nUnmappedTpl_);

      std::string tok = SafeFileToken(tgt.name);
      std::string thisMap = mapoutprefix.empty() ? std::string()
                                                 : mapoutprefix + tok + ".map";
      if (!mapoutprefix.empty()) {
        if (WriteMapFile(thisMap, *tgt.top, *parent.top, R, tgt.name, parent.name))
          return CpptrajState::ERR;
        mprintf("\t  Map written to '%s'\n", thisMap.c_str());
      }

      // TI dual files for every non-parent member (or all if tioutprefix set).
      if (!tioutprefix.empty() && (int)i != pIdx) {
        std::string tip = tioutprefix + tok;
        Frame tgtX, tplX;
        LoadCoords(*tgt.top, tgt.ds, tgtX);
        LoadCoords(*parent.top, parent.ds, tplX);
        Topology top0, top1;
        Frame frm0, frm1;
        NameType rn0 = ResNameOf(*parent.top, "L0");
        NameType rn1 = ResNameOf(*tgt.top, "L1");
        if (BuildTiUnit(true,  *tgt.top, *parent.top, tgtX, tplX, R.dual_, top0, frm0, rn0) ||
            BuildTiUnit(false, *tgt.top, *parent.top, tgtX, tplX, R.dual_, top1, frm1, rn1))
        {
          mprinterr("Error: timap: failed to build TI units for '%s'.\n", tgt.name.c_str());
          return CpptrajState::ERR;
        }
        std::string f0m = tip + ".0.mol2";
        std::string f1m = tip + ".1.mol2";
        std::string f0l = tip + ".0.lib";
        std::string f1l = tip + ".1.lib";
        std::string fsc = tip + ".scmask";
        std::string fat = tip + ".atoms";
        std::string u0 = OffUnitName(rn0.Truncated(), '0');
        std::string u1 = OffUnitName(rn1.Truncated(), '1');
        if (u0 == u1) { u0 += "0"; u1 += "1"; }
        top0.SetParmName(u0, FileName(f0m));
        top1.SetParmName(u1, FileName(f1m));
        if (WriteMol2File(f0m, top0, frm0, State.DSL())) return CpptrajState::ERR;
        if (WriteMol2File(f1m, top1, frm1, State.DSL())) return CpptrajState::ERR;
        if (WriteAmberLib(f0l, u0, top0, frm0)) return CpptrajState::ERR;
        if (WriteAmberLib(f1l, u1, top1, frm1)) return CpptrajState::ERR;
        if (WriteScmask(fsc, top0, top1, R.dual_)) return CpptrajState::ERR;
        if (WriteDualAtoms(fat, top0, top1, R.dual_)) return CpptrajState::ERR;
        if (mapoutprefix.empty()) {
          std::string fmap = tip + ".map";
          if (WriteMapFile(fmap, *tgt.top, *parent.top, R, tgt.name, parent.name))
            return CpptrajState::ERR;
        }
        mprintf("\t  TI dual-topology prefix: %s\n", tip.c_str());
      }

      if (!writeLib) continue;
      Topology* aligned = tgt.top->ModifyByMap(R.outputOrder_);
      if (aligned == 0) {
        mprinterr("Error: timap: failed to apply atom order for '%s'.\n", tgt.name.c_str());
        return CpptrajState::ERR;
      }
      std::string thisOut = outprefix + tok + AlignedOutExt(outFmt);
      Frame tgtX, alignedX;
      LoadCoords(*tgt.top, tgt.ds, tgtX);
      ApplyOrderToFrame(tgtX, R.outputOrder_, alignedX);
      if (WriteAlignedOut(thisOut, outFmt, *aligned, alignedX, State.DSL())) {
        delete aligned;
        return CpptrajState::ERR;
      }
      mprintf("\t  Aligned %s: %s\n", AlignedOutLabel(outFmt), thisOut.c_str());
      delete aligned;
    }
    return CpptrajState::OK;
  }

  // ----- Single-pair / self-map path -----
  if (names.empty()) {
    mprinterr("Error: timap: no target topology specified.\n");
    return CpptrajState::ERR;
  }
  std::string tgtName = names[0];
  if (tplName.empty() && names.size() > 1)
    tplName = names[1];
  if (names.size() > 2) {
    mprinterr("Error: timap: extra arguments; use 'series' for more than one analog.\n");
    return CpptrajState::ERR;
  }

  DataSet* tgtDs = 0;
  Topology* tgt = ResolveTop(State, tgtName, &tgtDs);
  if (tgt == 0) {
    mprinterr("Error: timap: target '%s' not found.\n", tgtName.c_str());
    return CpptrajState::ERR;
  }
  Topology* tpl = tgt;
  DataSet* tplDs = tgtDs;
  std::string tplUsed = tgtName;
  if (!tplName.empty()) {
    tpl = ResolveTop(State, tplName, &tplDs);
    if (tpl == 0) {
      mprinterr("Error: timap: template '%s' not found.\n", tplName.c_str());
      return CpptrajState::ERR;
    }
    tplUsed = tplName;
  } else {
    tplName = tgtName;
    tplUsed = tgtName;
  }

  TIMatch matcher;
  matcher.SetDebug(State.Debug());
  matcher.SetUseNaOrder(naorder);
  matcher.SetUseAaOrder(aaorder);
  if (aaorder && anchor.empty())
    matcher.SetAnchorName("C");
  if (!anchor.empty()) matcher.SetAnchorName(anchor);
  if (!seedStr.empty()) {
    std::string l = ToLower(seedStr);
    if (l != "auto" && l != "names" && l != "na" && l != "aa" && l != "none") {
      mprinterr("Error: timap: unrecognized seed '%s'\n", seedStr.c_str());
      return CpptrajState::ERR;
    }
    matcher.SetSeed(TIMatch::SeedFromString(seedStr));
  }

  bool writeLib = !libout.empty() || !maponly;
  if (writeLib && libout.empty()) {
    libout = OffUnitName(ResNameOf(*tgt, "TGT").Truncated(), 'T');
    libout += AlignedOutExt(outFmt);
  }

  const char* walkMsg = "";
  if (naorder) walkMsg = " with nucleic-acid canonical walk";
  else if (aaorder) walkMsg = " with amino-acid canonical walk";
  bool selfMap = (tgt == tpl);
  if (selfMap)
    mprintf("    TIMAP: Reordering '%s' (%i atoms)%s.\n",
            tgtName.c_str(), tgt->Natom(), walkMsg);
  else
    mprintf("    TIMAP: Aligning '%s' (%i atoms) to template '%s' (%i atoms).\n",
            tgtName.c_str(), tgt->Natom(), tplUsed.c_str(), tpl->Natom());
  mprintf("\tSeed: %s\n", TIMatch::SeedStr(
            seedStr.empty() ? TIMatch::SEED_AUTO
                            : TIMatch::SeedFromString(seedStr)));
  if (naorder)
    mprintf("\tOrder: nucleic-acid canonical walk%s.\n",
            selfMap ? " of this residue" : " of the template");
  else if (aaorder)
    mprintf("\tOrder: amino-acid canonical walk (ff19SB)%s.\n",
            selfMap ? " of this residue" : " of the template");
  else if (selfMap)
    mprintf("\tOrder: this residue's current atom order (no template, no walk).\n");
  else
    mprintf("\tOrder: template file atom order.\n");
  mprintf("\tLeftover insertion anchor: %s\n",
          !anchor.empty() ? anchor.c_str() : (aaorder ? "C" : "O3'"));
  if (!tiout.empty())
    mprintf("\tTI dual-topology prefix: %s\n", tiout.c_str());
  if (writeLib)
    mprintf("\tAligned %s: %s\n", AlignedOutLabel(outFmt), libout.c_str());
  else if (maponly)
    mprintf("\tmaponly: not writing an aligned structure file.\n");

  TIMatch::Result R;
  if (matcher.Match(*tgt, *tpl, R)) return CpptrajState::ERR;

  mprintf("\tMapped %i / %i target atoms (%i insertions, %i template atoms unmatched).\n",
          R.nMapped_, tgt->Natom(), R.nInsertion_, R.nUnmappedTpl_);
  mprintf("\tTarget classified as '%s'; template classified as '%s'.\n",
          R.tgtKind_.c_str(), R.tplKind_.c_str());
  if (tgt->Natom() > 0 && R.outputOrder_.size() > 0) {
    int first = R.outputOrder_.front();
    int last  = R.outputOrder_.back();
    mprintf("\tAligned order starts at %s, ends at %s.\n",
            (*tgt)[first].c_str(), (*tgt)[last].c_str());
  }
  if (R.nMapped_ == 0)
    mprintf("Warning: No atoms were mapped. Check that template and target are related.\n");

  if (!mapout.empty()) {
    mprintf("\tMap written to '%s'\n", mapout.c_str());
    if (WriteMapFile(mapout, *tgt, *tpl, R, tgtName, tplUsed))
      return CpptrajState::ERR;
  }

  if (!tiout.empty()) {
    if (R.dual_.empty()) {
      mprinterr("Error: timap: empty dual-topology layout.\n");
      return CpptrajState::ERR;
    }
    Frame tgtX, tplX;
    bool haveTgtX = LoadCoords(*tgt, tgtDs, tgtX);
    bool haveTplX = LoadCoords(*tpl, tplDs, tplX);
    if (!haveTgtX || !haveTplX)
      mprintf("Warning: Missing coordinates for %s%s%s; dummy/real positions may be zero.\n",
              haveTgtX ? "" : "target",
              (!haveTgtX && !haveTplX) ? " and " : "",
              haveTplX ? "" : "template");

    Topology top0, top1;
    Frame frm0, frm1;
    NameType rn0 = ResNameOf(*tpl, "L0");
    NameType rn1 = ResNameOf(*tgt, "L1");
    if (BuildTiUnit(true,  *tgt, *tpl, tgtX, tplX, R.dual_, top0, frm0, rn0) ||
        BuildTiUnit(false, *tgt, *tpl, tgtX, tplX, R.dual_, top1, frm1, rn1))
    {
      mprinterr("Error: timap: failed to build TI units.\n");
      return CpptrajState::ERR;
    }
    if (top0.Natom() != top1.Natom()) {
      mprinterr("Error: timap: lambda-0/1 atom counts differ (%i vs %i).\n",
                top0.Natom(), top1.Natom());
      return CpptrajState::ERR;
    }
    mprintf("\tDual topology: %i atoms (%i shared, %i dummy in lambda 0, %i dummy in lambda 1).\n",
            top0.Natom(), R.nMapped_, R.nInsertion_, R.nUnmappedTpl_);

    std::string f0m = tiout + ".0.mol2";
    std::string f1m = tiout + ".1.mol2";
    std::string f0l = tiout + ".0.lib";
    std::string f1l = tiout + ".1.lib";
    std::string fsc = tiout + ".scmask";
    std::string fat = tiout + ".atoms";
    std::string u0 = OffUnitName(rn0.Truncated(), '0');
    std::string u1 = OffUnitName(rn1.Truncated(), '1');
    if (u0 == u1) {
      u0 += "0";
      u1 += "1";
    }
    top0.SetParmName(u0, FileName(f0m));
    top1.SetParmName(u1, FileName(f1m));
    if (WriteMol2File(f0m, top0, frm0, State.DSL())) return CpptrajState::ERR;
    if (WriteMol2File(f1m, top1, frm1, State.DSL())) return CpptrajState::ERR;
    if (WriteAmberLib(f0l, u0, top0, frm0)) return CpptrajState::ERR;
    if (WriteAmberLib(f1l, u1, top1, frm1)) return CpptrajState::ERR;
    if (WriteScmask(fsc, top0, top1, R.dual_)) return CpptrajState::ERR;
    if (WriteDualAtoms(fat, top0, top1, R.dual_)) return CpptrajState::ERR;
    if (mapout.empty()) {
      std::string fmap = tiout + ".map";
      mprintf("\tMap written to '%s'\n", fmap.c_str());
      if (WriteMapFile(fmap, *tgt, *tpl, R, tgtName, tplUsed))
        return CpptrajState::ERR;
    }
  }

  if (!writeLib && newname.empty() && !replace)
    return CpptrajState::OK;

  Topology* aligned = tgt->ModifyByMap(R.outputOrder_);
  if (aligned == 0) {
    mprinterr("Error: timap: failed to apply atom order.\n");
    return CpptrajState::ERR;
  }
  if (writeLib) {
    Frame tgtX, alignedX;
    LoadCoords(*tgt, tgtDs, tgtX);
    ApplyOrderToFrame(tgtX, R.outputOrder_, alignedX);
    if (WriteAlignedOut(libout, outFmt, *aligned, alignedX, State.DSL())) {
      delete aligned;
      return CpptrajState::ERR;
    }
  }

  if (replace) {
    if (tgtDs != 0 && tgtDs->Type() == DataSet::TOPOLOGY) {
      ((DataSet_Topology*)tgtDs)->SetTop(*aligned);
      mprintf("\tReplaced topology '%s' with aligned atom order.\n", tgtName.c_str());
    } else {
      mprinterr("Error: timap: replace requires <tgt> to be a topology set.\n");
      delete aligned;
      return CpptrajState::ERR;
    }
  }
  if (!newname.empty()) {
    aligned->SetParmName(newname, tgt->OriginalFilename());
    if (State.AddTopology(*aligned, newname)) {
      delete aligned;
      return CpptrajState::ERR;
    }
    mprintf("\tAligned topology added as '%s'\n", newname.c_str());
  }
  delete aligned;
  return CpptrajState::OK;
}
