#include "Exec_TIMatch.h"
#include "TemplateMatch.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "DataSet_Topology.h"
#include "DataSet_Coords.h"
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
#include <vector>

/* Source encoding: UTF-8. Unicode (O5′, 2′, λ) is comments-only.
 *
 * Workflow:
 *   1. Match target onto template (partial map).
 *   2. Write the target as an Amber OFF .lib in template atom order (out).
 *   3. If tiout <prefix> is given, also build a dual-topology pair of the
 *      same size: unique atoms are real in one end state and charge-0 /
 *      mass-0 / type-DUM dummies in the other.
 */

// Exec_TIMatch::Help()
void Exec_TIMatch::Help() const {
  mprintf("\t<tgt> [template <name>] [out <file>] [tiout <prefix>]\n"
          "\t[mapout <file>] [name <newparm>] [maponly] [replace] [naorder]\n"
          "\t[seed {auto|names|na|none}] [anchor <atomname>]\n"
          "  Align atoms in topology <tgt> to a user-supplied template so that\n"
          "  shared atoms occupy the same indices. Intended for thermodynamic\n"
          "  integration (TI) of nucleotides, amino acids, and small molecules.\n"
          "  Partial maps are expected (unlike atommap).\n"
          "  The template is lambda=0; <tgt> is lambda=1.\n"
          "  Default: write <tgt> as an Amber OFF .lib in template atom order\n"
          "  (out <file>; default <residue>.sorted.lib). No dummy atoms.\n"
          "  'maponly' skips that library (and in-memory parm) but still writes\n"
          "  mapout / tiout files if those keywords are given.\n"
          "  'tiout <prefix>' additionally writes a dual-topology TI pair:\n"
          "    <prefix>.0.mol2 / <prefix>.0.lib  lambda=0 template + dummy insertions\n"
          "    <prefix>.1.mol2 / <prefix>.1.lib  lambda=1 target + dummy unmatched atoms\n"
          "    <prefix>.scmask                   pmemd scmask1 / scmask2\n"
          "    <prefix>.atoms                    per-slot kind, names, and charges\n"
          "  Dummy atoms copy the partner's name and coords, have charge 0, mass 0,\n"
          "  and Amber type DUM.\n"
          "  Both TI files have the same atom count so residue indices match in pmemd.\n"
          "  'name' stores a remapped topology in memory; 'replace' overwrites <tgt>.\n"
          "  If 'template' is omitted, <tgt> is matched to itself (use with 'naorder'\n"
          "  to freeze a nucleic-acid template from an existing residue).\n"
          "  Official ModXNA parent fragments ship in $CPPTRAJHOME/dat/templatematch/.\n"
          "  Aliases: templatematch, timap.\n"
          "\n"
          "  Complete run (Amber OFF libraries parent.lib and analog.lib).\n"
          "  readdata, parm, and timatch run when entered (immediate commands).\n"
          "  Use go if the input also has trajin/actions; it is safe to include either way:\n"
          "    > readdata parent.lib name parent\n"
          "    > readdata analog.lib name analog\n"
          "    > timatch analog[analog] template parent[parent] out analog.lib\n"
          "    > go\n"
          "  Use readdata (not parm) for .lib files. The COORDS set is Name[Unit];\n"
          "  if the unit inside parent.lib is not 'parent', use parent[UnitName].\n"
          "  For nucleic acids add naorder so the shared-atom walk is\n"
          "  P -> OP -> O5' -> C5' -> C4' -> O4' -> C1' -> base -> C3' -> C2' -> O3':\n"
          "    > timatch analog[analog] template parent[parent] naorder out analog.lib\n"
          "    > go\n"
          "  Dual-topology TI (opt-in; dummy atoms, matching NATOM):\n"
          "    > timatch analog[analog] template parent[parent] naorder tiout analog_ti\n"
          "    > go\n"
          "  Mol2 inputs instead of OFF:\n"
          "    > parm parent.mol2 name parent\n"
          "    > parm analog.mol2 name analog\n"
          "    > timatch analog template parent out analog.lib\n"
          "    > go\n");
}

/** Look up a topology by the name the user typed.
  *
  * Order of attempts:
  *   1. A TOPOLOGY dataset (parm / parmwrite sets).
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
                       std::vector<TemplateMatch::Result::DualSlot> const& dual,
                       Topology& outTop, Frame& outFrm,
                       NameType const& resname)
{
  typedef TemplateMatch::Result::DualSlot Slot;
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
        mprintf("Warning: timatch: dummy %s has no bonded neighbor in the dual layout.\n",
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
        mprintf("Warning: timatch: dummy %s has no bonded neighbor in the dual layout.\n",
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
    mprinterr("Error: timatch: could not set up mol2 '%s'\n", fname.c_str());
    return 1;
  }
  if (out.WriteSingle(0, frm)) {
    mprinterr("Error: timatch: writing mol2 '%s'\n", fname.c_str());
    return 1;
  }
  out.EndTraj();
  mprintf("\tWrote '%s' (%i atoms).\n", fname.c_str(), top.Natom());
  return 0;
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

static int WriteScmask(std::string const& fname,
                       Topology const& top0, Topology const& top1,
                       std::vector<TemplateMatch::Result::DualSlot> const& dual)
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
  out.Printf("# timatch dual-topology TI masks\n");
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
                          std::vector<TemplateMatch::Result::DualSlot> const& dual)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  out.Printf("# timatch dual-topology atoms\n");
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

/** Write the human-readable correspondence used by Test_TemplateMatch. */
static int WriteMapFile(std::string const& fname, Topology const& tgt, Topology const& tpl,
                        TemplateMatch::Result const& R, std::string const& tgtName,
                        std::string const& tplName)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open map file '%s'\n", fname.c_str());
    return 1;
  }
  out.Printf("# timatch tgt='%s' template='%s'\n", tgtName.c_str(), tplName.c_str());
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

/** Parse arguments, run TemplateMatch::Match, write an aligned OFF library,
  * and optionally dual-topology TI files.
  */
// Exec_TIMatch::Execute()
Exec::RetType Exec_TIMatch::Execute(CpptrajState& State, ArgList& argIn) {
  std::string mapout = argIn.GetStringKey("mapout");
  std::string libout = argIn.GetStringKey("out");
  std::string newname = argIn.GetStringKey("name");
  std::string tplName = argIn.GetStringKey("template");
  std::string seedStr = argIn.GetStringKey("seed");
  std::string anchor = argIn.GetStringKey("anchor");
  std::string tiout = argIn.GetStringKey("tiout");
  bool maponly = argIn.hasKey("maponly");
  bool replace = argIn.hasKey("replace");
  bool naorder = argIn.hasKey("naorder");

  std::string tgtName = argIn.GetStringNext();
  if (tgtName.empty()) {
    mprinterr("Error: timatch: no target topology specified.\n");
    return CpptrajState::ERR;
  }
  if (tplName.empty())
    tplName = argIn.GetStringNext();

  DataSet* tgtDs = 0;
  Topology* tgt = FindNamedTop(State.DSL(), tgtName, &tgtDs);
  if (tgt == 0) {
    mprinterr("Error: timatch: target '%s' not found.\n", tgtName.c_str());
    return CpptrajState::ERR;
  }
  Topology* tpl = tgt;
  DataSet* tplDs = tgtDs;
  std::string tplUsed = tgtName;
  if (!tplName.empty()) {
    tpl = FindNamedTop(State.DSL(), tplName, &tplDs);
    if (tpl == 0) {
      mprinterr("Error: timatch: template '%s' not found.\n", tplName.c_str());
      return CpptrajState::ERR;
    }
    tplUsed = tplName;
  } else {
    tplName = tgtName;
    tplUsed = tgtName;
  }

  TemplateMatch matcher;
  matcher.SetDebug(State.Debug());
  matcher.SetUseNaOrder(naorder);
  if (!anchor.empty()) matcher.SetAnchorName(anchor);
  if (!seedStr.empty()) {
    std::string l = ToLower(seedStr);
    if (l != "auto" && l != "names" && l != "na" && l != "none") {
      mprinterr("Error: timatch: unrecognized seed '%s'\n", seedStr.c_str());
      return CpptrajState::ERR;
    }
    matcher.SetSeed(TemplateMatch::SeedFromString(seedStr));
  }

  bool writeLib = !libout.empty() || !maponly;
  if (writeLib && libout.empty()) {
    libout = OffUnitName(ResNameOf(*tgt, "TGT").Truncated(), 'T');
    libout += ".sorted.lib";
  }

  mprintf("    TIMATCH: Aligning '%s' (%i atoms) to template '%s' (%i atoms).\n",
          tgtName.c_str(), tgt->Natom(), tplUsed.c_str(), tpl->Natom());
  mprintf("\tSeed: %s\n", TemplateMatch::SeedStr(
            seedStr.empty() ? TemplateMatch::SEED_AUTO
                            : TemplateMatch::SeedFromString(seedStr)));
  if (naorder)
    mprintf("\tUsing nucleic-acid canonical walk as the template order.\n");
  else
    mprintf("\tUsing template file atom order as the shared-atom order.\n");
  mprintf("\tLeftover insertion anchor: %s\n",
          anchor.empty() ? "O3'" : anchor.c_str());
  if (!tiout.empty())
    mprintf("\tTI dual-topology prefix: %s\n", tiout.c_str());
  if (writeLib)
    mprintf("\tAligned OFF library: %s\n", libout.c_str());
  else if (maponly)
    mprintf("\tmaponly: not writing an aligned OFF library.\n");

  TemplateMatch::Result R;
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
      mprinterr("Error: timatch: empty dual-topology layout.\n");
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
      mprinterr("Error: timatch: failed to build TI units.\n");
      return CpptrajState::ERR;
    }
    if (top0.Natom() != top1.Natom()) {
      mprinterr("Error: timatch: lambda-0/1 atom counts differ (%i vs %i).\n",
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
    mprinterr("Error: timatch: failed to apply atom order.\n");
    return CpptrajState::ERR;
  }
  if (writeLib) {
    Frame tgtX, alignedX;
    LoadCoords(*tgt, tgtDs, tgtX);
    ApplyOrderToFrame(tgtX, R.outputOrder_, alignedX);
    std::string unit = OffUnitName(ResNameOf(*aligned, "TGT").Truncated(), 'T');
    aligned->SetParmName(unit, FileName(libout));
    if (WriteAmberLib(libout, unit, *aligned, alignedX)) {
      delete aligned;
      return CpptrajState::ERR;
    }
  }

  if (replace) {
    if (tgtDs != 0 && tgtDs->Type() == DataSet::TOPOLOGY) {
      ((DataSet_Topology*)tgtDs)->SetTop(*aligned);
      mprintf("\tReplaced topology '%s' with aligned atom order.\n", tgtName.c_str());
    } else {
      mprinterr("Error: timatch: replace requires <tgt> to be a topology set.\n");
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
