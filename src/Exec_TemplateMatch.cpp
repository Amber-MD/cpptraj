#include "Exec_TemplateMatch.h"
#include "TemplateMatch.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "DataSet_Topology.h"
#include "DataSet_Coords.h"
#include "StringRoutines.h"
#include "ArgList.h"
#include <algorithm>
#include <vector>

/* Source encoding: UTF-8. Unicode (O5′, 2′, λ) is comments-only; Help() text
 * and string literals stay ASCII. This Exec looks up topologies, runs
 * TemplateMatch, writes mapout, and optionally stores a remapped parm.
 */

void Exec_TemplateMatch::Help() const {
  mprintf("\t<tgt> [template <name>] [mapout <file>] [name <newparm>]\n"
          "\t[maponly] [replace] [naorder]\n"
          "\t[seed {auto|names|na|none}] [anchor <atomname>]\n"
          "  Align atoms in topology <tgt> to a user-supplied <template> so that\n"
          "  shared atoms occupy the same indices. Intended for thermodynamic\n"
          "  integration (TI) of nucleotides, amino acids, and small molecules.\n"
          "  Partial maps are expected (unlike atommap).\n"
          "  The template atom order *is* the TI shared-atom order. Extra atoms on\n"
          "  <tgt> are inserted next to the atom they are bonded to. Completely\n"
          "  unmatched leftovers go just before <anchor> (default O3').\n"
          "  If 'template' is omitted, <tgt> is matched to itself (use with 'naorder'\n"
          "  to freeze a nucleic-acid template from an existing residue).\n"
          "  'naorder' walks the template as P -> O5' -> C5' -> C4' -> O4' -> C1' ->\n"
          "  base -> C3' -> C2' -> O3' instead of using the file atom order.\n"
          "  'seed auto' uses nucleic-acid scaffold roles when they are detected,\n"
          "  otherwise unique atom names.\n"
          "  Official ModXNA parent fragments ship in $CPPTRAJHOME/dat/templatematch/.\n"
          "  Generating your own template (pseudocode):\n"
          "    1. Pick the lambda=0 parent (unmodified dA, benzene, cysteine, ...).\n"
          "    2. Load it. For nucleotides, optionally:\n"
          "         templatematch <parent> naorder name <parent>.ti\n"
          "       For ligands / amino acids the file order is the TI order.\n"
          "    3. Map each analog:\n"
          "         templatematch <analog> template <parent> mapout analog.map name analog.ti\n"
          "    4. Amber OFF (.lib via readdata) is required when element cannot be\n"
          "       guessed from the atom name (e.g. SE is sulfur in mol2; use elmnt 34).\n");
}

/** Look up a topology by the name the user typed.
  *
  * Order of attempts:
  *   1. A TOPOLOGY dataset (parm / parmwrite sets).
  *   2. A COORDS dataset — Amber OFF units loaded with readdata appear as
  *      LibName[UnitName], e.g. FLE[FLE] or CYS[CYS].
  *   3. A purely numeric string, treated as a parm index.
  *
  * \param dsOut If non-null, receives the TOPOLOGY dataset pointer when the
  *              hit was a topology set (needed for the 'replace' keyword).
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

/** Write the human-readable correspondence used by Test_TemplateMatch.
  *
  * Columns: output index, original target atom, template partner (or ---),
  * insertion flag. Unmapped template atoms are listed after the table so
  * a partial map (2′-OH vs H2″, OH vs H1, Se vs S) is obvious in the file.
  */
static int WriteMapFile(std::string const& fname, Topology const& tgt, Topology const& tpl,
                        TemplateMatch::Result const& R, std::string const& tgtName,
                        std::string const& tplName)
{
  CpptrajFile out;
  if (out.OpenWrite(fname)) {
    mprinterr("Error: Could not open map file '%s'\n", fname.c_str());
    return 1;
  }
  out.Printf("# templatematch tgt='%s' template='%s'\n", tgtName.c_str(), tplName.c_str());
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

/** Parse arguments, run TemplateMatch::Match, optionally write the map and
  * store a remapped topology. 'maponly' skips the topology write so tests
  * can check the correspondence without mutating parm sets.
  */
Exec::RetType Exec_TemplateMatch::Execute(CpptrajState& State, ArgList& argIn) {
  std::string mapout = argIn.GetStringKey("mapout");
  std::string newname = argIn.GetStringKey("name");
  std::string tplName = argIn.GetStringKey("template");
  std::string seedStr = argIn.GetStringKey("seed");
  std::string anchor = argIn.GetStringKey("anchor");
  bool maponly = argIn.hasKey("maponly");
  bool replace = argIn.hasKey("replace");
  bool naorder = argIn.hasKey("naorder");

  std::string tgtName = argIn.GetStringNext();
  if (tgtName.empty()) {
    mprinterr("Error: templatematch: no target topology specified.\n");
    return CpptrajState::ERR;
  }
  if (tplName.empty())
    tplName = argIn.GetStringNext();

  DataSet* tgtDs = 0;
  Topology* tgt = FindNamedTop(State.DSL(), tgtName, &tgtDs);
  if (tgt == 0) {
    mprinterr("Error: templatematch: target '%s' not found.\n", tgtName.c_str());
    return CpptrajState::ERR;
  }
  Topology* tpl = tgt;
  std::string tplUsed = tgtName;
  if (!tplName.empty()) {
    tpl = FindNamedTop(State.DSL(), tplName, 0);
    if (tpl == 0) {
      mprinterr("Error: templatematch: template '%s' not found.\n", tplName.c_str());
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
      mprinterr("Error: templatematch: unrecognized seed '%s'\n", seedStr.c_str());
      return CpptrajState::ERR;
    }
    matcher.SetSeed(TemplateMatch::SeedFromString(seedStr));
  }

  mprintf("    TEMPLATEMATCH: Aligning '%s' (%i atoms) to template '%s' (%i atoms).\n",
          tgtName.c_str(), tgt->Natom(), tplUsed.c_str(), tpl->Natom());
  mprintf("\tSeed: %s\n", TemplateMatch::SeedStr(
            seedStr.empty() ? TemplateMatch::SEED_AUTO
                            : TemplateMatch::SeedFromString(seedStr)));
  if (naorder)
    mprintf("\tUsing nucleic-acid canonical walk as the template order.\n");
  else
    mprintf("\tUsing template file atom order as the TI shared-atom order.\n");
  mprintf("\tLeftover insertion anchor: %s\n",
          anchor.empty() ? "O3'" : anchor.c_str());

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

  if (maponly)
    return CpptrajState::OK;

  Topology* aligned = tgt->ModifyByMap(R.outputOrder_);
  if (aligned == 0) {
    mprinterr("Error: templatematch: failed to apply atom order.\n");
    return CpptrajState::ERR;
  }
  aligned->SetParmName(newname.empty() ? (tgtName + ".aligned") : newname,
                       tgt->OriginalFilename());

  if (replace) {
    if (tgtDs != 0 && tgtDs->Type() == DataSet::TOPOLOGY) {
      ((DataSet_Topology*)tgtDs)->SetTop(*aligned);
      mprintf("\tReplaced topology '%s' with aligned atom order.\n", tgtName.c_str());
    } else {
      mprinterr("Error: replace requires <tgt> to be a topology set.\n");
      delete aligned;
      return CpptrajState::ERR;
    }
  }
  if (!newname.empty()) {
    if (State.AddTopology(*aligned, newname)) {
      delete aligned;
      return CpptrajState::ERR;
    }
    mprintf("\tAligned topology added as '%s'\n", newname.c_str());
  } else if (!replace) {
    std::string autoName = tgtName + ".aligned";
    if (State.AddTopology(*aligned, autoName)) {
      delete aligned;
      return CpptrajState::ERR;
    }
    mprintf("\tAligned topology added as '%s'\n", autoName.c_str());
  }
  delete aligned;
  return CpptrajState::OK;
}
