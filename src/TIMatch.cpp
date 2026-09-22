#include "TemplateMatch.h"
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

/* Source encoding: UTF-8.
 *
 * Unicode in this file is confined to comments (O5′, C1′, 2′, χ, λ). Every
 * string literal and identifier stays ASCII so MSVC compiles the file without
 * /utf-8. Do not put a UTF-8 prime inside "quotes" that the compiler sees.
 *
 * TemplateMatch implements a partial, graph-based correspondence from a
 * user target residue onto a user template residue. The template's current
 * atom order is the thermodynamic-integration (λ = 0 / λ = 1 shared) order
 * unless the caller requested naorder, in which case the template is walked
 * as a nucleotide:
 *
 *   P → non-bridging OP → O5′ → C5′ → C4′ → O4′ → C1′ →
 *   nucleobase starting at the glycosidic nitrogen χ →
 *   C3′ → C2′ and its 2′ substituents → O3′
 *
 * Three chemistries are handled by the same code:
 *   - Nucleic acids: scaffold roles seed the map (FLE rA vs ERN dA).
 *   - Amino acids: unique names + grow (Cys vs selenocysteine; S ≠ Se).
 *   - Small molecules: unique names + grow (benzene vs phenol).
 *
 * Nested Graph / Scaffold types live in an anonymous namespace so the free
 * helpers below can use them without making them TemplateMatch members.
 */

// -----------------------------------------------------------------------------
TemplateMatch::TemplateMatch() :
  seed_(SEED_AUTO),
  debug_(0),
  useNaOrder_(false),
  useAaOrder_(false),
  anchorName_("O3'")
{}

const char* TemplateMatch::SeedStr(SeedType t) {
  switch (t) {
    case SEED_AUTO:  return "auto";
    case SEED_NAMES: return "names";
    case SEED_NA:    return "na";
    case SEED_AA:    return "aa";
    case SEED_NONE:  return "none";
  }
  return "auto";
}

/// Parse seed {auto|names|na|none}; unknown strings fall back to auto.
TemplateMatch::SeedType TemplateMatch::SeedFromString(std::string const& s) {
  std::string l = ToLower(s);
  if (l == "auto")  return SEED_AUTO;
  if (l == "names") return SEED_NAMES;
  if (l == "na")    return SEED_NA;
  if (l == "aa")    return SEED_AA;
  if (l == "none")  return SEED_NONE;
  return SEED_AUTO;
}

namespace {

typedef std::vector<int> Iarray;

// -----------------------------------------------------------------------------
/// True for elemental hydrogen (not extra points / Drude particles).
static inline bool IsH(Atom const& a) {
  return (a.Element() == Atom::HYDROGEN);
}

/// Heavy = not H, extra point, or Drude. Used for ring detection and signatures.
static inline bool IsHeavy(Atom const& a) {
  return (!IsH(a) && a.Element() != Atom::EXTRAPT && a.Element() != Atom::DRUDE);
}

/// Truncated Amber name (no padding). Comparisons and map files use this.
static std::string AtomName(Atom const& a) {
  return a.Name().Truncated();
}

// -----------------------------------------------------------------------------
/// Bonded graph over one Topology: full, heavy-only, and hydrogen neighbor lists.
/** Built once per residue so DetectScaffold / Grow / CanonicalWalk can ask
  * "who is bonded to C2′?" without walking Atom::bond_iterator each time.
  */
class Graph {
  public:
    Graph() : top_(0) {}
    explicit Graph(Topology const& t) { Setup(t); }

    /// Fill nbr_ / heavyNbr_ / hNbr_ from Topology bonds. Indices match atom #.
    void Setup(Topology const& t) {
      top_ = &t;
      int n = t.Natom();
      nbr_.assign(n, Iarray());
      heavyNbr_.assign(n, Iarray());
      hNbr_.assign(n, Iarray());
      for (int i = 0; i < n; i++) {
        Atom const& at = t[i];
        for (Atom::bond_iterator b = at.bondbegin(); b != at.bondend(); ++b) {
          nbr_[i].push_back(*b);
          if (IsHeavy(t[*b]))
            heavyNbr_[i].push_back(*b);
          else if (IsH(t[*b]))
            hNbr_[i].push_back(*b);
        }
      }
    }

    Topology const& Top() const { return *top_; }
    int Natom()            const { return top_->Natom(); }
    Atom const& operator[](int i) const { return (*top_)[i]; }
    std::string Name(int i) const { return AtomName((*top_)[i]); }
    Atom::AtomicElementType Elt(int i) const { return (*top_)[i].Element(); }
    bool IsHydrogen(int i) const { return IsH((*top_)[i]); }
    int HeavyDeg(int i)    const { return (int)heavyNbr_[i].size(); }
    Iarray const& Nbr(int i)      const { return nbr_[i]; }
    Iarray const& HeavyNbr(int i) const { return heavyNbr_[i]; }
    Iarray const& HNbr(int i)     const { return hNbr_[i]; }

    /// Every atom whose truncated name equals \p name (0, 1, or many hits).
    Iarray FindByName(std::string const& name) const {
      Iarray hits;
      for (int i = 0; i < Natom(); i++) {
        if (Name(i) == name)
          hits.push_back(i);
      }
      return hits;
    }

    Iarray HeavyNbrsOf(int i) const { return heavyNbr_[i]; }

  private:
    Topology const* top_;
    std::vector<Iarray> nbr_;
    std::vector<Iarray> heavyNbr_;
    std::vector<Iarray> hNbr_;
};

// -----------------------------------------------------------------------------
/// Nucleic-acid (or ligand) roles detected on one residue.
/** role_[P] … role_[CHI] are atom indices or −1. sugarRing_ is the chosen
  * furanose (4 C + 1 O/S). base_ is the connected component hanging off χ,
  * excluding C1′. kind_ is nucleotide | sugar | base | unknown.
  *
  * Amino acids fill aa_ (N, CA, C, O, CB) and kind_ "amino" when the
  * peptide N–CA–C=O motif is found. Small molecules stay kind "unknown".
  */
class Scaffold {
  public:
    /// Canonical nucleotide roles. CHI is the glycosidic N (or C) of χ.
    enum Role { P = 0, O5p, C5p, C4p, O4p, C1p, C2p, C3p, O3p, CHI, NROLES };
    enum AaRole { AA_N = 0, AA_CA, AA_C, AA_O, AA_CB, AA_NROLES };
    static const char* RoleStr(Role r) {
      static const char* n[NROLES] = {
        "P", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O3'", "chi"
      };
      return n[r];
    }
    static const char* AaRoleStr(AaRole r) {
      static const char* n[AA_NROLES] = { "N", "CA", "C", "O", "CB" };
      return n[r];
    }

    Scaffold() : hasP_(false), ok_(false), aaOk_(false), isPro_(false) {
      for (int i = 0; i < NROLES; i++) role_[i] = -1;
      for (int i = 0; i < AA_NROLES; i++) aa_[i] = -1;
      kind_ = "unknown";
      family_ = "none";
      twoPrime_ = "unknown";
    }

    int Get(Role r) const { return role_[r]; }
    void Set(Role r, int i) { role_[r] = i; }
    int Aa(AaRole r) const { return aa_[r]; }
    void SetAa(AaRole r, int i) { aa_[r] = i; }
    bool HasSugar() const { return !sugarRing_.empty(); }
    /// Combined nucleotide: furanose + glycosidic N + a nucleobase component.
    bool IsCombined() const {
      return HasSugar() && role_[CHI] >= 0 && !base_.empty();
    }

    int role_[NROLES];
    Iarray op_;
    Iarray sugarRing_;
    Iarray base_;
    std::string kind_;
    std::string family_;
    std::string twoPrime_;
    bool hasP_;
    bool ok_;
    bool aaOk_;
    bool isPro_;
    int aa_[AA_NROLES];
    std::vector<std::string> notes_;
};

// -----------------------------------------------------------------------------
typedef std::vector<Iarray> Cycles;

/// Rotate a cycle so the lowest index is first; pick the lexicographically
/// smaller of the two directions so the same ring is not stored twice.
static Iarray NormalizeCycle(Iarray path) {
  if (path.empty()) return path;
  int i0 = 0;
  for (int i = 1; i < (int)path.size(); i++) {
    if (path[i] < path[i0]) i0 = i;
  }
  Iarray rot;
  rot.reserve(path.size());
  for (int k = 0; k < (int)path.size(); k++)
    rot.push_back(path[(i0 + k) % path.size()]);
  if (rot.size() > 2 && rot.back() < rot[1]) {
    Iarray rev(1, rot[0]);
    for (int k = (int)rot.size() - 1; k >= 1; k--)
      rev.push_back(rot[k]);
    int j0 = 0;
    for (int i = 1; i < (int)rev.size(); i++) {
      if (rev[i] < rev[j0]) j0 = i;
    }
    Iarray rot2;
    for (int k = 0; k < (int)rev.size(); k++)
      rot2.push_back(rev[(j0 + k) % rev.size()]);
    return rot2;
  }
  return rot;
}

/// Depth-first search for simple heavy-atom cycles of a fixed length.
/** \p start is the origin of this search (only neighbors ≥ start are
  * followed, which together with NormalizeCycle deduplicates rings).
  */
static void DfsCycle(int start, Iarray const& path, std::set<int> const& seen,
                     std::vector<Iarray> const& heavyNbr, int length,
                     std::set<Iarray>& found)
{
  if ((int)path.size() == length) {
    int last = path.back();
    Iarray const& nbr = heavyNbr[last];
    if (std::find(nbr.begin(), nbr.end(), start) != nbr.end() && path[0] == start)
      found.insert(NormalizeCycle(path));
    return;
  }
  int last = path.back();
  Iarray const& nbr = heavyNbr[last];
  for (Iarray::const_iterator it = nbr.begin(); it != nbr.end(); ++it) {
    int nxt = *it;
    if (seen.find(nxt) != seen.end()) continue;
    if (nxt < start) continue;
    Iarray path2 = path;
    path2.push_back(nxt);
    std::set<int> seen2 = seen;
    seen2.insert(nxt);
    DfsCycle(start, path2, seen2, heavyNbr, length, found);
  }
}

/// All unique heavy-atom cycles of \p length (5 = furanose / imidazole, 6 = pyrimidine).
static Cycles CyclesOfLength(Graph const& g, int length) {
  std::set<Iarray> found;
  std::vector<Iarray> heavy;
  heavy.reserve(g.Natom());
  for (int i = 0; i < g.Natom(); i++)
    heavy.push_back(g.HeavyNbr(i));
  for (int i = 0; i < g.Natom(); i++) {
    Iarray path(1, i);
    std::set<int> seen;
    seen.insert(i);
    DfsCycle(i, path, seen, heavy, length, found);
  }
  Cycles out(found.begin(), found.end());
  return out;
}

/// Count O / C / N / S among the atoms of one cycle (used to score sugar vs nucleobase rings).
static void CountRingElts(Graph const& g, Iarray const& cyc,
                          int& nO, int& nC, int& nN, int& nS)
{
  nO = nC = nN = nS = 0;
  for (Iarray::const_iterator it = cyc.begin(); it != cyc.end(); ++it) {
    switch (g.Elt(*it)) {
      case Atom::OXYGEN:   nO++; break;
      case Atom::CARBON:   nC++; break;
      case Atom::NITROGEN: nN++; break;
      case Atom::SULFUR:   nS++; break;
      default: break;
    }
  }
}

/// Furanose: four carbons and exactly one O or S, no nitrogen (so not imidazole).
static bool IsSugarRing(Graph const& g, Iarray const& cyc) {
  int nO, nC, nN, nS;
  CountRingElts(g, cyc, nO, nC, nN, nS);
  return (nC == 4 && nN == 0 && (nO + nS) == 1);
}

/// Prefer the 5-ring whose carbons have the most nucleotide-like exo substituents
/// (glycosidic N, C5′ carbon, O3′/O5′ oxygen). Ties break lexicographically.
static Iarray PickSugarRing(Graph const& g, Cycles const& rings5) {
  int bestScore = -1;
  Iarray best;
  for (Cycles::const_iterator c = rings5.begin(); c != rings5.end(); ++c) {
    if (!IsSugarRing(g, *c)) continue;
    std::set<int> ring(c->begin(), c->end());
    int score = 0;
    for (Iarray::const_iterator it = c->begin(); it != c->end(); ++it) {
      if (g.Elt(*it) != Atom::CARBON) continue;
      bool exoC = false, exoO = false, exoN = false;
      Iarray const& hn = g.HeavyNbr(*it);
      for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
        if (ring.find(*j) != ring.end()) continue;
        if (g.Elt(*j) == Atom::CARBON)   exoC = true;
        if (g.Elt(*j) == Atom::OXYGEN)   exoO = true;
        if (g.Elt(*j) == Atom::NITROGEN) exoN = true;
      }
      if (exoC) score += 3;
      if (exoO) score += 2;
      if (exoN) score += 4;
    }
    if (score > bestScore || (score == bestScore && (best.empty() || *c < best))) {
      bestScore = score;
      best = *c;
    }
  }
  return best;
}

/// Connected component of χ, stopping at C1′ so the sugar is not swallowed into the base.
static Iarray BaseAtoms(Graph const& g, int chi, int c1) {
  std::set<int> seen;
  seen.insert(chi);
  if (c1 >= 0) seen.insert(c1);
  Iarray stack(1, chi);
  Iarray out;
  while (!stack.empty()) {
    int i = stack.back();
    stack.pop_back();
    out.push_back(i);
    Iarray const& nbr = g.Nbr(i);
    for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
      if (seen.find(*j) != seen.end()) continue;
      if (c1 >= 0 && *j == c1) continue;
      seen.insert(*j);
      stack.push_back(*j);
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

/// purine if an imidazole 5-ring and a 6-ring with N are both present; else pyrimidine.
static std::string BaseFamily(Graph const& g, int chi,
                              Cycles const& rings5, Cycles const& rings6)
{
  bool imidazole = false;
  for (Cycles::const_iterator c = rings5.begin(); c != rings5.end(); ++c) {
    int nO, nC, nN, nS;
    CountRingElts(g, *c, nO, nC, nN, nS);
    if (nN >= 1 && nO == 0) imidazole = true;
  }
  bool sixN = false;
  for (Cycles::const_iterator c = rings6.begin(); c != rings6.end(); ++c) {
    int nO, nC, nN, nS;
    CountRingElts(g, *c, nO, nC, nN, nS);
    if (nN >= 1) sixN = true;
  }
  if (imidazole && sixN) return "purine";
  if (sixN) return "pyrimidine";
  if (imidazole) return "purine";
  int nC = 0;
  Iarray const& hn = g.HeavyNbr(chi);
  for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
    if (g.Elt(*j) == Atom::CARBON) nC++;
  }
  if (nC >= 2) return "purine";
  return "pyrimidine";
}

/// Classify the 2′ substituent: HH (deoxy), HOH (ribose), OMe, OR, F, other.
static std::string TwoPrimePattern(Graph const& g,
                                   Scaffold const& sc)
{
  int c2 = sc.Get(Scaffold::C2p);
  if (c2 < 0) return "unknown";
  std::set<int> skip;
  if (sc.Get(Scaffold::C1p) >= 0) skip.insert(sc.Get(Scaffold::C1p));
  if (sc.Get(Scaffold::C3p) >= 0) skip.insert(sc.Get(Scaffold::C3p));
  if (sc.Get(Scaffold::O4p) >= 0) skip.insert(sc.Get(Scaffold::O4p));
  int nO = 0, nC = 0, nF = 0, nH = (int)g.HNbr(c2).size();
  int oIdx = -1;
  Iarray const& hn = g.HeavyNbr(c2);
  for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
    if (skip.find(*j) != skip.end()) continue;
    if (g.Elt(*j) == Atom::OXYGEN) { nO++; oIdx = *j; }
    else if (g.Elt(*j) == Atom::CARBON) nC++;
    else if (g.Elt(*j) == Atom::FLUORINE) nF++;
  }
  if (nO == 0 && nC == 0 && nF == 0 && nH >= 1) return "HH";
  if (nO == 1 && nC == 0 && oIdx >= 0) {
    if ((int)g.HeavyNbr(oIdx).size() == 1) return "HOH";
    Iarray carbons;
    Iarray const& ohn = g.HeavyNbr(oIdx);
    for (Iarray::const_iterator j = ohn.begin(); j != ohn.end(); ++j) {
      if (*j != c2 && g.Elt(*j) == Atom::CARBON)
        carbons.push_back(*j);
    }
    if (carbons.size() == 1 && g.HeavyDeg(carbons[0]) == 1) return "OMe";
    if (!carbons.empty()) return "OR";
    return "OR";
  }
  if (nF == 1) return "F";
  return "other";
}

/// Best glycosidic nitrogen for a nucleobase fragment (no sugar): N9 / N1 names win.
static int GuessGlycosidicN(Graph const& g,
                            Cycles const& rings5, Cycles const& rings6)
{
  int best = -1;
  int bestScore = -999;
  for (int i = 0; i < g.Natom(); i++) {
    if (g.Elt(i) != Atom::NITROGEN) continue;
    int nC = 0;
    Iarray const& hn = g.HeavyNbr(i);
    for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
      if (g.Elt(*j) == Atom::CARBON) nC++;
    }
    bool in5 = false, in6 = false;
    for (Cycles::const_iterator c = rings5.begin(); c != rings5.end(); ++c) {
      if (std::find(c->begin(), c->end(), i) == c->end()) continue;
      int nO, nC2, nN, nS;
      CountRingElts(g, *c, nO, nC2, nN, nS);
      if (nN >= 1 && nO == 0) in5 = true;
    }
    for (Cycles::const_iterator c = rings6.begin(); c != rings6.end(); ++c) {
      if (std::find(c->begin(), c->end(), i) != c->end()) in6 = true;
    }
    int score = nC;
    if (in5) score += 8;
    if (in6 && !in5) score += 3;
    std::string nm = ToLower(g.Name(i));
    if (nm == "n9") score += 10;
    else if (nm == "n1" && !in5) score += 10;
    if (score > bestScore) {
      bestScore = score;
      best = i;
    }
  }
  return best;
}

/// Fill any still-empty roles from unique Amber names (P, O5′, C1′, N9, …).
static void NameHints(Graph const& g, Scaffold& sc) {
  static const char* names[] = {
    "P", "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O3'", "N9", "N1", "S4'", 0
  };
  static const Scaffold::Role roles[] = {
    Scaffold::P, Scaffold::O5p, Scaffold::C5p,
    Scaffold::C4p, Scaffold::O4p, Scaffold::C1p,
    Scaffold::C2p, Scaffold::C3p, Scaffold::O3p,
    Scaffold::CHI, Scaffold::CHI, Scaffold::O4p
  };
  for (int k = 0; names[k] != 0; k++) {
    Scaffold::Role role = roles[k];
    if (sc.Get(role) >= 0) continue;
    Iarray hits = g.FindByName(names[k]);
    if (hits.size() == 1)
      sc.Set(role, hits[0]);
  }
}

/// If the topology carries ModXNA HEAD/TAIL/ANCHOR metadata, use it as χ / O3′ / C1′.
static void ApplyModxna(Graph const& g, Scaffold& sc) {
  ModXNA_Info const& mx = g.Top().Modxna();
  if (!mx.HasModxna()) return;
  std::string head = mx.Head();
  if (!head.empty() && head[0] == '@') head = head.substr(1);
  std::string tail = mx.Tail();
  if (!tail.empty() && tail[0] == '@') tail = tail.substr(1);
  std::string anchor = mx.Anchor();
  if (!anchor.empty() && anchor[0] == '@') anchor = anchor.substr(1);
  Iarray hh = head.empty() ? Iarray() : g.FindByName(head);
  Iarray tt = tail.empty() ? Iarray() : g.FindByName(tail);
  Iarray aa = anchor.empty() ? Iarray() : g.FindByName(anchor);
  int hi = (hh.size() == 1) ? hh[0] : -1;
  int ti = (tt.size() == 1) ? tt[0] : -1;
  int ai = (aa.size() == 1) ? aa[0] : -1;
  if (hi >= 0 && sc.Get(Scaffold::CHI) < 0 &&
      (g.Elt(hi) == Atom::NITROGEN || g.Elt(hi) == Atom::CARBON))
    sc.Set(Scaffold::CHI, hi);
  if (ti >= 0 && sc.Get(Scaffold::O3p) < 0 && g.Elt(ti) == Atom::OXYGEN)
    sc.Set(Scaffold::O3p, ti);
  if (ai >= 0 && sc.Get(Scaffold::C1p) < 0)
    sc.Set(Scaffold::C1p, ai);
}

/// Assign C1′/C4′ (the two carbons on O4′), then C2′/C3′, C5′, O5′, OP, O3′.
/** C1′ is the O4′-bonded carbon that also bears the glycosidic nitrogen (or
  * the less O-rich exo carbon). C4′ is the other O4′-bonded carbon and should
  * lead to C5′/O5′.
  */
static void AssignSugarRoles(Graph const& g, Scaffold& sc,
                             Iarray const& sugar)
{
  std::set<int> ring(sugar.begin(), sugar.end());
  int o4 = -1;
  Iarray ringC;
  for (Iarray::const_iterator it = sugar.begin(); it != sugar.end(); ++it) {
    if (g.Elt(*it) == Atom::OXYGEN || g.Elt(*it) == Atom::SULFUR) o4 = *it;
    if (g.Elt(*it) == Atom::CARBON) ringC.push_back(*it);
  }
  sc.Set(Scaffold::O4p, o4);
  Iarray bondedToO;
  if (o4 >= 0) {
    for (Iarray::const_iterator it = ringC.begin(); it != ringC.end(); ++it) {
      Iarray const& hn = g.HeavyNbr(*it);
      if (std::find(hn.begin(), hn.end(), o4) != hn.end())
        bondedToO.push_back(*it);
    }
  }
  int c1 = -1, c4 = -1;
  if (bondedToO.size() == 2) {
    std::vector< std::pair<int,int> > scored;
    for (Iarray::const_iterator it = bondedToO.begin(); it != bondedToO.end(); ++it) {
      int score = 0;
      Iarray const& hn = g.HeavyNbr(*it);
      for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
        if (ring.find(*j) != ring.end()) continue;
        if (g.Elt(*j) == Atom::NITROGEN) score -= 10;
        else if (g.Elt(*j) == Atom::OXYGEN) score += 2;
        else if (g.Elt(*j) == Atom::CARBON) {
          int nO = 0;
          Iarray const& jn = g.HeavyNbr(*j);
          for (Iarray::const_iterator k = jn.begin(); k != jn.end(); ++k)
            if (g.Elt(*k) == Atom::OXYGEN) nO++;
          if (nO) score += 6;
          else if ((int)jn.size() <= 1) score -= 4;
          else score += 1;
        }
      }
      scored.push_back(std::make_pair(score, *it));
    }
    std::sort(scored.begin(), scored.end());
    c1 = scored[0].second;
    c4 = scored[1].second;
  } else if (bondedToO.size() == 1) {
    c1 = bondedToO[0];
  }
  sc.Set(Scaffold::C1p, c1);
  sc.Set(Scaffold::C4p, c4);
  if (c1 < 0 || c4 < 0) {
    sc.notes_.push_back("could not split C1'/C4' on the sugar ring");
    return;
  }
  Iarray remain;
  for (Iarray::const_iterator it = ringC.begin(); it != ringC.end(); ++it) {
    if (*it != c1 && *it != c4) remain.push_back(*it);
  }
  int c3 = -1;
  for (Iarray::const_iterator it = remain.begin(); it != remain.end(); ++it) {
    Iarray const& hn = g.HeavyNbr(*it);
    if (std::find(hn.begin(), hn.end(), c4) != hn.end()) { c3 = *it; break; }
  }
  int c2 = -1;
  for (Iarray::const_iterator it = remain.begin(); it != remain.end(); ++it) {
    if (*it != c3) { c2 = *it; break; }
  }
  if (c2 >= 0 && c1 >= 0) {
    Iarray const& hn = g.HeavyNbr(c2);
    bool c2c1 = std::find(hn.begin(), hn.end(), c1) != hn.end();
    if (!c2c1 && c3 >= 0) {
      Iarray const& hn3 = g.HeavyNbr(c3);
      if (std::find(hn3.begin(), hn3.end(), c1) != hn3.end())
        std::swap(c2, c3);
    }
  }
  sc.Set(Scaffold::C3p, c3);
  sc.Set(Scaffold::C2p, c2);

  if (c4 >= 0) {
    Iarray exoC;
    Iarray const& hn = g.HeavyNbr(c4);
    for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
      if (ring.find(*j) == ring.end() && g.Elt(*j) == Atom::CARBON)
        exoC.push_back(*j);
    }
    if (!exoC.empty()) {
      int c5 = exoC[0];
      sc.Set(Scaffold::C5p, c5);
      Iarray o5cand;
      Iarray const& c5n = g.HeavyNbr(c5);
      for (Iarray::const_iterator j = c5n.begin(); j != c5n.end(); ++j) {
        if (g.Elt(*j) == Atom::OXYGEN) o5cand.push_back(*j);
      }
      int p = sc.Get(Scaffold::P);
      if (p >= 0) {
        int o5 = -1;
        for (Iarray::const_iterator j = o5cand.begin(); j != o5cand.end(); ++j) {
          Iarray const& jn = g.HeavyNbr(*j);
          if (std::find(jn.begin(), jn.end(), p) != jn.end()) { o5 = *j; break; }
        }
        sc.Set(Scaffold::O5p, o5 >= 0 ? o5 : (o5cand.empty() ? -1 : o5cand[0]));
      } else if (!o5cand.empty()) {
        sc.Set(Scaffold::O5p, o5cand[0]);
      }
    }
  }

  int p = sc.Get(Scaffold::P);
  int o5 = sc.Get(Scaffold::O5p);
  if (p >= 0) {
    Iarray ops;
    Iarray const& pn = g.HeavyNbr(p);
    for (Iarray::const_iterator j = pn.begin(); j != pn.end(); ++j) {
      if (g.Elt(*j) == Atom::OXYGEN && *j != o5) ops.push_back(*j);
    }
    std::sort(ops.begin(), ops.end(), [&g](int a, int b){ return g.Name(a) < g.Name(b); });
    sc.op_ = ops;
  }

  if (c3 >= 0) {
    Iarray o3cand;
    Iarray const& hn = g.HeavyNbr(c3);
    for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
      if (ring.find(*j) == ring.end() && g.Elt(*j) == Atom::OXYGEN)
        o3cand.push_back(*j);
    }
    if (!o3cand.empty())
      sc.Set(Scaffold::O3p, o3cand[0]);
  }
}

/// Unique truncated name, or −1 if missing / not unique.
static int UniqueName(Graph const& g, const char* nm) {
  Iarray hits = g.FindByName(nm);
  return (hits.size() == 1) ? hits[0] : -1;
}

/// Terminal oxygen (heavy degree 1): carbonyl O, carboxylate O, hydroxyl O.
static bool IsTermO(Graph const& g, int i) {
  return (g.Elt(i) == Atom::OXYGEN && g.HeavyDeg(i) == 1);
}

/// Fill peptide N / CA / C / O / CB from Amber names, then the N–CA–C=O motif.
/** ff19SB amino19.lib is the baseline: N, H, CA, HA, side chain from CB, C, O
  * (proline: N, CD … CB, CA, C, O). Noncanonical residues that keep backbone
  * names are labeled the same way; graph backup handles renamed backbones.
  */
static void DetectAmino(Graph const& g, Scaffold& sc) {
  sc.SetAa(Scaffold::AA_N,  UniqueName(g, "N"));
  sc.SetAa(Scaffold::AA_CA, UniqueName(g, "CA"));
  sc.SetAa(Scaffold::AA_C,  UniqueName(g, "C"));
  sc.SetAa(Scaffold::AA_O,  UniqueName(g, "O"));
  sc.SetAa(Scaffold::AA_CB, UniqueName(g, "CB"));

  if (sc.Aa(Scaffold::AA_N) < 0 || sc.Aa(Scaffold::AA_CA) < 0 ||
      sc.Aa(Scaffold::AA_C) < 0) {
    for (int ni = 0; ni < g.Natom(); ni++) {
      if (g.Elt(ni) != Atom::NITROGEN) continue;
      Iarray const& nNbr = g.HeavyNbr(ni);
      for (Iarray::const_iterator ca = nNbr.begin(); ca != nNbr.end(); ++ca) {
        if (g.Elt(*ca) != Atom::CARBON) continue;
        Iarray const& caNbr = g.HeavyNbr(*ca);
        for (Iarray::const_iterator cj = caNbr.begin(); cj != caNbr.end(); ++cj) {
          if (*cj == ni || g.Elt(*cj) != Atom::CARBON) continue;
          int oIdx = -1;
          int nTermO = 0;
          Iarray const& cNbr = g.HeavyNbr(*cj);
          for (Iarray::const_iterator oj = cNbr.begin(); oj != cNbr.end(); ++oj) {
            if (IsTermO(g, *oj)) { nTermO++; oIdx = *oj; }
          }
          if (nTermO != 1 || (int)cNbr.size() > 3) continue;
          if (sc.Aa(Scaffold::AA_N) < 0)  sc.SetAa(Scaffold::AA_N, ni);
          if (sc.Aa(Scaffold::AA_CA) < 0) sc.SetAa(Scaffold::AA_CA, *ca);
          if (sc.Aa(Scaffold::AA_C) < 0)  sc.SetAa(Scaffold::AA_C, *cj);
          if (sc.Aa(Scaffold::AA_O) < 0)  sc.SetAa(Scaffold::AA_O, oIdx);
        }
      }
    }
  }

  int n  = sc.Aa(Scaffold::AA_N);
  int ca = sc.Aa(Scaffold::AA_CA);
  int c  = sc.Aa(Scaffold::AA_C);
  if (sc.Aa(Scaffold::AA_O) < 0 && c >= 0) {
    Iarray const& cn = g.HeavyNbr(c);
    for (Iarray::const_iterator j = cn.begin(); j != cn.end(); ++j) {
      if (IsTermO(g, *j)) { sc.SetAa(Scaffold::AA_O, *j); break; }
    }
  }
  if (sc.Aa(Scaffold::AA_CB) < 0 && ca >= 0) {
    Iarray const& can = g.HeavyNbr(ca);
    for (Iarray::const_iterator j = can.begin(); j != can.end(); ++j) {
      if (g.Elt(*j) == Atom::CARBON && *j != c)
        sc.SetAa(Scaffold::AA_CB, *j);
    }
  }

  n  = sc.Aa(Scaffold::AA_N);
  ca = sc.Aa(Scaffold::AA_CA);
  c  = sc.Aa(Scaffold::AA_C);
  sc.aaOk_ = (n >= 0 && ca >= 0 && c >= 0);
  if (!sc.aaOk_) return;

  sc.isPro_ = false;
  if (n >= 0) {
    int nC = 0;
    Iarray const& nn = g.HeavyNbr(n);
    for (Iarray::const_iterator j = nn.begin(); j != nn.end(); ++j) {
      if (g.Elt(*j) == Atom::CARBON) nC++;
    }
    sc.isPro_ = (nC >= 2);
  }

  if (sc.kind_ == "unknown") {
    sc.kind_ = "amino";
    if (sc.isPro_) sc.family_ = "pro";
    else if (sc.Aa(Scaffold::AA_CB) < 0) sc.family_ = "gly";
    else sc.family_ = "std";
  }
}

/// Detect nucleotide / sugar / base / amino / unknown and fill Scaffold roles.
/** Amino acids without a furanose get kind "amino" (N–CA–C=O). Ligands
  * (benzene, phenol) stay "unknown" and unique-name seed.
  */
void DetectScaffold(Graph const& g, Scaffold& sc) {
  sc = Scaffold();
  if (g.Natom() == 0) {
    sc.notes_.push_back("empty topology");
    return;
  }
  Cycles rings5 = CyclesOfLength(g, 5);
  Cycles rings6 = CyclesOfLength(g, 6);

  Iarray pAtoms;
  for (int i = 0; i < g.Natom(); i++) {
    if (g.Elt(i) == Atom::PHOSPHORUS) pAtoms.push_back(i);
  }
  if (!pAtoms.empty()) {
    sc.Set(Scaffold::P, pAtoms[0]);
    sc.hasP_ = true;
    if (pAtoms.size() > 1)
      sc.notes_.push_back("multiple P atoms; using first");
  }

  Iarray sugar = PickSugarRing(g, rings5);
  if (!sugar.empty()) {
    sc.sugarRing_ = sugar;
    AssignSugarRoles(g, sc, sugar);
  } else {
    sc.notes_.push_back("no furanose ring");
  }

  int c1 = sc.Get(Scaffold::C1p);
  if (c1 >= 0) {
    Iarray nchi;
    Iarray const& hn = g.HeavyNbr(c1);
    for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
      if (g.Elt(*j) == Atom::NITROGEN) nchi.push_back(*j);
    }
    if (!nchi.empty()) {
      sc.Set(Scaffold::CHI, nchi[0]);
      sc.base_ = BaseAtoms(g, nchi[0], c1);
      sc.family_ = BaseFamily(g, nchi[0], rings5, rings6);
    }
  }

  if (sc.Get(Scaffold::C2p) >= 0)
    sc.twoPrime_ = TwoPrimePattern(g, sc);

  if (sc.Get(Scaffold::CHI) >= 0 && !sugar.empty()) {
    sc.kind_ = "nucleotide";
    sc.ok_ = true;
  } else if (!sugar.empty()) {
    sc.kind_ = "sugar";
    sc.ok_ = (sc.Get(Scaffold::C1p) >= 0 && sc.Get(Scaffold::O4p) >= 0);
  } else {
    int nN = 0, nC = 0;
    for (int i = 0; i < g.Natom(); i++) {
      if (g.Elt(i) == Atom::NITROGEN) nN++;
      if (g.Elt(i) == Atom::CARBON) nC++;
    }
    bool looksBase = (nN >= 2 && nC >= 4 && (!rings5.empty() || !rings6.empty()));
    if (sc.Get(Scaffold::CHI) >= 0 || looksBase) {
      sc.kind_ = "base";
      if (sc.Get(Scaffold::CHI) < 0) {
        int chi = GuessGlycosidicN(g, rings5, rings6);
        sc.Set(Scaffold::CHI, chi);
        if (chi >= 0) {
          sc.base_ = BaseAtoms(g, chi, -1);
          sc.family_ = BaseFamily(g, chi, rings5, rings6);
        }
      }
      sc.ok_ = (sc.Get(Scaffold::CHI) >= 0);
    } else {
      sc.kind_ = "unknown";
    }
  }

  NameHints(g, sc);
  ApplyModxna(g, sc);

  int chi = sc.Get(Scaffold::CHI);
  if (chi >= 0) {
    sc.base_ = BaseAtoms(g, chi, sc.Get(Scaffold::C1p));
    sc.family_ = BaseFamily(g, chi, rings5, rings6);
  }
  if (sc.Get(Scaffold::C2p) >= 0 && sc.kind_ != "base")
    sc.twoPrime_ = TwoPrimePattern(g, sc);

  // A C1′ methyl/base cap must not be labeled C5′ (common on nucleoside fragments).
  int c5 = sc.Get(Scaffold::C5p);
  c1 = sc.Get(Scaffold::C1p);
  int c4 = sc.Get(Scaffold::C4p);
  if (c5 >= 0 && c1 >= 0) {
    Iarray const& c1n = g.HeavyNbr(c1);
    bool onC1 = std::find(c1n.begin(), c1n.end(), c5) != c1n.end();
    bool onC4 = false;
    if (c4 >= 0) {
      Iarray const& c4n = g.HeavyNbr(c4);
      onC4 = std::find(c4n.begin(), c4n.end(), c5) != c4n.end();
    }
    if (onC1 && !onC4) {
      sc.Set(Scaffold::C5p, -1);
      sc.notes_.push_back("cleared C5' attached to C1'");
      c5 = -1;
    }
  }

  if (sc.kind_ == "unknown")
    DetectAmino(g, sc);
}

// -----------------------------------------------------------------------------
/// Neighbor signature used to grow the map: (element, heavy degree, sorted neighbor elements).
/** S vs Se (Cys/Sec) and O vs H (phenol/benzene at C1) fail this test, which
  * is how those substitutions become insertions rather than forced matches.
  */
struct AtomSig {
  Atom::AtomicElementType elt_;
  int deg_;
  std::vector<Atom::AtomicElementType> nbr_;
  bool operator<(AtomSig const& rhs) const {
    if (elt_ != rhs.elt_) return elt_ < rhs.elt_;
    if (deg_ != rhs.deg_) return deg_ < rhs.deg_;
    return nbr_ < rhs.nbr_;
  }
};

/// Build the grow signature of atom i (hydrogens are excluded from degree / neighbors).
static AtomSig MakeSig(Graph const& g, int i) {
  AtomSig s;
  s.elt_ = g.Elt(i);
  s.deg_ = g.HeavyDeg(i);
  Iarray const& hn = g.HeavyNbr(i);
  s.nbr_.reserve(hn.size());
  for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j)
    s.nbr_.push_back(g.Elt(*j));
  std::sort(s.nbr_.begin(), s.nbr_.end());
  return s;
}

/// Stable order for same-signature twins: unique names if possible, else atom index.
static Iarray OrderTwins(Graph const& g, Iarray kids) {
  std::set<std::string> names;
  for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
    names.insert(g.Name(*j));
  if ((int)names.size() == (int)kids.size()) {
    std::sort(kids.begin(), kids.end(), [&g](int a, int b){ return g.Name(a) < g.Name(b); });
    return kids;
  }
  std::sort(kids.begin(), kids.end());
  return kids;
}

/// Record tgt[n] ↔ tpl[r] unless either side is already claimed. Returns true if newly added.
static bool AddMap(Iarray& mapping, Iarray& usedTpl, int n, int r) {
  if (n < 0 || r < 0) return false;
  if (mapping[n] >= 0) return (mapping[n] == r);
  if (usedTpl[r] != 0) return false;
  mapping[n] = r;
  usedTpl[r] = 1;
  return true;
}

/// Grow the correspondence from already-mapped atoms by unique neighbor signatures.
/** Repeated until a pass adds nothing. Unique (1:1) signatures are paired;
  * same-size twin groups (H5′/H5″, OP1/OP2) are paired by name. This is the
  * step that places a 2′-OH next to C2′ without touching template H2″.
  */
void Grow(Graph const& tgt, Graph const& tpl,
          Iarray& mapping, Iarray& usedTpl)
{
  bool changed = true;
  while (changed) {
    changed = false;
    Iarray mappedTgt;
    for (int n = 0; n < tgt.Natom(); n++) {
      if (mapping[n] >= 0) mappedTgt.push_back(n);
    }
    for (Iarray::const_iterator it = mappedTgt.begin(); it != mappedTgt.end(); ++it) {
      int n = *it;
      int r = mapping[n];
      Iarray nUn, rUn;
      Iarray const& nn = tgt.Nbr(n);
      for (Iarray::const_iterator j = nn.begin(); j != nn.end(); ++j) {
        if (mapping[*j] < 0) nUn.push_back(*j);
      }
      Iarray const& rn = tpl.Nbr(r);
      for (Iarray::const_iterator j = rn.begin(); j != rn.end(); ++j) {
        if (!usedTpl[*j]) rUn.push_back(*j);
      }
      if (nUn.empty() || rUn.empty()) continue;
      std::map<AtomSig, Iarray> nBy, rBy;
      for (Iarray::const_iterator j = nUn.begin(); j != nUn.end(); ++j)
        nBy[MakeSig(tgt, *j)].push_back(*j);
      for (Iarray::const_iterator j = rUn.begin(); j != rUn.end(); ++j)
        rBy[MakeSig(tpl, *j)].push_back(*j);
      for (std::map<AtomSig, Iarray>::iterator nb = nBy.begin(); nb != nBy.end(); ++nb) {
        std::map<AtomSig, Iarray>::iterator rb = rBy.find(nb->first);
        if (rb == rBy.end()) continue;
        Iarray& nlist = nb->second;
        Iarray& rlist = rb->second;
        if (nlist.size() == 1 && rlist.size() == 1) {
          if (AddMap(mapping, usedTpl, nlist[0], rlist[0]))
            changed = true;
        } else if (nlist.size() == rlist.size() && nlist.size() > 1) {
          Iarray ns = OrderTwins(tgt, nlist);
          Iarray rs = OrderTwins(tpl, rlist);
          for (size_t k = 0; k < ns.size(); k++) {
            if (AddMap(mapping, usedTpl, ns[k], rs[k]))
              changed = true;
          }
        }
      }
    }
  }
}

/// Pair leftover atoms that share a unique name *and* the same element.
/** Element check is required: mol2 atom name SE would otherwise look like
  * sulfur, and even with OFF, SE must not be forced onto SG.
  */
void MatchUniqueNames(Graph const& tgt, Graph const& tpl,
                      Iarray& mapping, Iarray& usedTpl)
{
  std::map<std::string, Iarray> refByName;
  for (int i = 0; i < tpl.Natom(); i++)
    refByName[tpl.Name(i)].push_back(i);
  for (int i = 0; i < tgt.Natom(); i++) {
    if (mapping[i] >= 0) continue;
    std::map<std::string, Iarray>::const_iterator it = refByName.find(tgt.Name(i));
    if (it == refByName.end()) continue;
    Iarray hits;
    for (Iarray::const_iterator r = it->second.begin(); r != it->second.end(); ++r) {
      if (!usedTpl[*r]) hits.push_back(*r);
    }
    if (hits.size() == 1 && tgt.Elt(i) == tpl.Elt(hits[0]))
      AddMap(mapping, usedTpl, i, hits[0]);
  }
}

/// Zip hydrogens of already-mapped heavies, ordered by name (H61 before H62).
void MatchHydrogens(Graph const& tgt, Graph const& tpl,
                    Iarray& mapping, Iarray& usedTpl)
{
  for (int n = 0; n < tgt.Natom(); n++) {
    int r = mapping[n];
    if (r < 0) continue;
    Iarray nH, rH;
    Iarray const& th = tgt.HNbr(n);
    for (Iarray::const_iterator j = th.begin(); j != th.end(); ++j) {
      if (mapping[*j] < 0) nH.push_back(*j);
    }
    Iarray const& rh = tpl.HNbr(r);
    for (Iarray::const_iterator j = rh.begin(); j != rh.end(); ++j) {
      if (!usedTpl[*j]) rH.push_back(*j);
    }
    if (nH.empty() || rH.empty()) continue;
    Iarray ns = OrderTwins(tgt, nH);
    Iarray rs = OrderTwins(tpl, rH);
    size_t nzip = std::min(ns.size(), rs.size());
    for (size_t k = 0; k < nzip; k++)
      AddMap(mapping, usedTpl, ns[k], rs[k]);
  }
}

/// Seed, grow, unique-name, grow, hydrogens. mapping[tgt] = tpl or −1.
/** SEED_AUTO uses nucleic-acid roles when either residue looks like NA,
  * amino-acid roles when the peptide N–CA–C=O motif is found, otherwise names.
  */
void MapGraphs(Graph const& tgt, Graph const& tpl,
               Scaffold const& tgtSc, Scaffold const& tplSc,
               Iarray& mapping, TemplateMatch::SeedType seedIn)
{
  mapping.assign(tgt.Natom(), -1);
  Iarray usedTpl(tpl.Natom(), 0);

  TemplateMatch::SeedType mode = seedIn;
  if (mode == TemplateMatch::SEED_AUTO) {
    bool na = tgtSc.ok_ || tplSc.ok_ ||
              tgtSc.Get(Scaffold::P) >= 0 || tplSc.Get(Scaffold::P) >= 0 ||
              tgtSc.Get(Scaffold::C1p) >= 0 || tplSc.Get(Scaffold::C1p) >= 0 ||
              tgtSc.Get(Scaffold::CHI) >= 0 || tplSc.Get(Scaffold::CHI) >= 0;
    bool aa = tgtSc.aaOk_ || tplSc.aaOk_;
    if (na)      mode = TemplateMatch::SEED_NA;
    else if (aa) mode = TemplateMatch::SEED_AA;
    else         mode = TemplateMatch::SEED_NAMES;
  }

  if (mode == TemplateMatch::SEED_NA) {
    static const Scaffold::Role SEED[] = {
      Scaffold::P, Scaffold::O5p, Scaffold::C5p, Scaffold::C4p, Scaffold::O4p,
      Scaffold::C1p, Scaffold::C3p, Scaffold::C2p, Scaffold::O3p, Scaffold::CHI
    };
    for (int k = 0; k < 10; k++)
      AddMap(mapping, usedTpl, tgtSc.Get(SEED[k]), tplSc.Get(SEED[k]));
    size_t nop = std::min(tgtSc.op_.size(), tplSc.op_.size());
    for (size_t i = 0; i < nop; i++)
      AddMap(mapping, usedTpl, tgtSc.op_[i], tplSc.op_[i]);
  }

  if (mode == TemplateMatch::SEED_AA) {
    static const Scaffold::AaRole SEED[] = {
      Scaffold::AA_N, Scaffold::AA_CA, Scaffold::AA_C, Scaffold::AA_O, Scaffold::AA_CB
    };
    for (int k = 0; k < 5; k++)
      AddMap(mapping, usedTpl, tgtSc.Aa(SEED[k]), tplSc.Aa(SEED[k]));
  }

  if (mode != TemplateMatch::SEED_NONE) {
    Grow(tgt, tpl, mapping, usedTpl);
    MatchUniqueNames(tgt, tpl, mapping, usedTpl);
    Grow(tgt, tpl, mapping, usedTpl);
    MatchHydrogens(tgt, tpl, mapping, usedTpl);
  }
}

// -----------------------------------------------------------------------------
/// Depth-first emit of an unmapped insertion subtree (heavy atoms before hydrogens).
static void EmitTree(Graph const& g, int n, int parent,
                     std::vector<char>& used, Iarray const& mappedNew,
                     Iarray& order)
{
  if (used[n]) return;
  used[n] = 1;
  order.push_back(n);
  Iarray kids;
  Iarray const& nbr = g.Nbr(n);
  for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
    if (*j == parent || used[*j] || mappedNew[*j]) continue;
    kids.push_back(*j);
  }
  std::sort(kids.begin(), kids.end(), [&g](int a, int b) {
    if (g.IsHydrogen(a) != g.IsHydrogen(b)) return !g.IsHydrogen(a) && g.IsHydrogen(b);
    if (g.Elt(a) != g.Elt(b)) return g.Elt(a) < g.Elt(b);
    return g.Name(a) < g.Name(b);
  });
  for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
    EmitTree(g, *j, n, used, mappedNew, order);
}

/// Permutation of target atoms: walk the template, emit the mapped atom, then insertions.
/** Insertions bonded to a mapped parent (O2′ on C2′, OH on C1, Se on CB) are
  * placed immediately after that parent. Completely unmatched leftovers go
  * just before O3′ (or \p anchorName). Unmapped template atoms are skipped
  * (H2″ of dA, H1 of benzene, SG of cysteine).
  */
Iarray OutputOrder(Graph const& tgt, Graph const& tpl,
                   Scaffold const& tplSc,
                   Iarray const& mapping,
                   Iarray const& parentOrder,
                   std::string const& anchorName)
{
  Iarray refToNew(tpl.Natom(), -1);
  Iarray mappedNew(tgt.Natom(), 0);
  for (int n = 0; n < tgt.Natom(); n++) {
    if (mapping[n] >= 0) {
      refToNew[mapping[n]] = n;
      mappedNew[n] = 1;
    }
  }
  std::vector<char> used(tgt.Natom(), 0);
  Iarray order;
  int o3slot = -1;
  int o3 = tplSc.Get(Scaffold::O3p);
  if (o3 < 0 && !anchorName.empty()) {
    Iarray hits = tpl.FindByName(anchorName);
    if (hits.size() == 1) o3 = hits[0];
  }

  for (Iarray::const_iterator pr = parentOrder.begin(); pr != parentOrder.end(); ++pr) {
    int r = *pr;
    if (r == o3) o3slot = (int)order.size();
    if (r < 0 || r >= tpl.Natom() || refToNew[r] < 0) continue;
    int n = refToNew[r];
    if (used[n]) continue;
    used[n] = 1;
    order.push_back(n);
    Iarray kids;
    Iarray const& nbr = tgt.Nbr(n);
    for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
      if (!used[*j] && !mappedNew[*j]) kids.push_back(*j);
    }
    std::sort(kids.begin(), kids.end(), [&tgt](int a, int b) {
      if (tgt.IsHydrogen(a) != tgt.IsHydrogen(b))
        return !tgt.IsHydrogen(a) && tgt.IsHydrogen(b);
      if (tgt.Elt(a) != tgt.Elt(b)) return tgt.Elt(a) < tgt.Elt(b);
      return tgt.Name(a) < tgt.Name(b);
    });
    for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
      EmitTree(tgt, *j, n, used, mappedNew, order);
  }

  Iarray leftover;
  for (int i = 0; i < tgt.Natom(); i++) {
    if (!used[i]) leftover.push_back(i);
  }
  std::sort(leftover.begin(), leftover.end(), [&tgt](int a, int b) {
    if (tgt.IsHydrogen(a) != tgt.IsHydrogen(b))
      return !tgt.IsHydrogen(a) && tgt.IsHydrogen(b);
    if (tgt.Elt(a) != tgt.Elt(b)) return tgt.Elt(a) < tgt.Elt(b);
    if (tgt.Name(a) != tgt.Name(b)) return tgt.Name(a) < tgt.Name(b);
    return a < b;
  });
  if (!leftover.empty()) {
    int insertAt = (o3slot >= 0) ? o3slot : (int)order.size();
    order.insert(order.begin() + insertAt, leftover.begin(), leftover.end());
  }
  return order;
}

/// Dual-topology layout: every template atom keeps a slot; insertions are extra slots.
/** SHARED: both real. TPL_ONLY: real on λ=0, dummy on λ=1. TGT_ONLY: dummy on λ=0,
  * real on λ=1. Slot order matches OutputOrder for target atoms, with unmatched
  * template atoms inserted at their parent-order positions (not skipped).
  */
static std::vector<TemplateMatch::Result::DualSlot>
DualOrder(Graph const& tgt, Graph const& tpl,
          Scaffold const& tplSc,
          Iarray const& mapping,
          Iarray const& parentOrder,
          std::string const& anchorName)
{
  typedef TemplateMatch::Result::DualSlot Slot;
  Iarray refToNew(tpl.Natom(), -1);
  Iarray mappedNew(tgt.Natom(), 0);
  for (int n = 0; n < tgt.Natom(); n++) {
    if (mapping[n] >= 0) {
      refToNew[mapping[n]] = n;
      mappedNew[n] = 1;
    }
  }
  std::vector<char> used(tgt.Natom(), 0);
  std::vector<Slot> dual;
  int o3slot = -1;
  int o3 = tplSc.Get(Scaffold::O3p);
  if (o3 < 0 && !anchorName.empty()) {
    Iarray hits = tpl.FindByName(anchorName);
    if (hits.size() == 1) o3 = hits[0];
  }

  for (Iarray::const_iterator pr = parentOrder.begin(); pr != parentOrder.end(); ++pr) {
    int r = *pr;
    if (r == o3) o3slot = (int)dual.size();
    if (r < 0 || r >= tpl.Natom()) continue;
    if (refToNew[r] < 0) {
      dual.push_back(Slot(r, -1));
      continue;
    }
    int n = refToNew[r];
    if (used[n]) continue;
    used[n] = 1;
    dual.push_back(Slot(r, n));
    Iarray kids;
    Iarray const& nbr = tgt.Nbr(n);
    for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
      if (!used[*j] && !mappedNew[*j]) kids.push_back(*j);
    }
    std::sort(kids.begin(), kids.end(), [&tgt](int a, int b) {
      if (tgt.IsHydrogen(a) != tgt.IsHydrogen(b))
        return !tgt.IsHydrogen(a) && tgt.IsHydrogen(b);
      if (tgt.Elt(a) != tgt.Elt(b)) return tgt.Elt(a) < tgt.Elt(b);
      return tgt.Name(a) < tgt.Name(b);
    });
    Iarray tree;
    for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
      EmitTree(tgt, *j, n, used, mappedNew, tree);
    for (Iarray::const_iterator j = tree.begin(); j != tree.end(); ++j)
      dual.push_back(Slot(-1, *j));
  }

  Iarray leftover;
  for (int i = 0; i < tgt.Natom(); i++) {
    if (!used[i]) leftover.push_back(i);
  }
  std::sort(leftover.begin(), leftover.end(), [&tgt](int a, int b) {
    if (tgt.IsHydrogen(a) != tgt.IsHydrogen(b))
      return !tgt.IsHydrogen(a) && tgt.IsHydrogen(b);
    if (tgt.Elt(a) != tgt.Elt(b)) return tgt.Elt(a) < tgt.Elt(b);
    if (tgt.Name(a) != tgt.Name(b)) return tgt.Name(a) < tgt.Name(b);
    return a < b;
  });
  if (!leftover.empty()) {
    int insertAt = (o3slot >= 0) ? o3slot : (int)dual.size();
    std::vector<Slot> extra;
    extra.reserve(leftover.size());
    for (Iarray::const_iterator j = leftover.begin(); j != leftover.end(); ++j)
      extra.push_back(Slot(-1, *j));
    dual.insert(dual.begin() + insertAt, extra.begin(), extra.end());
  }
  return dual;
}

/// Mark atom i used and append it; no-op if i < 0 or already emitted.
static void EmitIdx(int i, std::vector<char>& used, Iarray& order) {
  if (i < 0 || used[i]) return;
  used[i] = 1;
  order.push_back(i);
}

/// Emit a heavy atom and its hydrogens (name-sorted) for the canonical NA walk.
static void EmitHeavyAndH(Graph const& g, int i,
                          std::vector<char>& used, Iarray& order)
{
  EmitIdx(i, used, order);
  if (i < 0) return;
  Iarray hs = g.HNbr(i);
  std::sort(hs.begin(), hs.end(), [&g](int a, int b){ return g.Name(a) < g.Name(b); });
  for (Iarray::const_iterator h = hs.begin(); h != hs.end(); ++h)
    EmitIdx(*h, used, order);
}

/// Lower number = earlier in the nucleobase BFS (ring atoms before exo O/N/C, H last).
static int BasePriority(Graph const& g, int i) {
  if (g.IsHydrogen(i)) return 9;
  int nRingish = g.HeavyDeg(i);
  if (g.Elt(i) == Atom::OXYGEN && nRingish <= 1) return 5;
  if (g.Elt(i) == Atom::NITROGEN && nRingish <= 1) return 4;
  if (g.Elt(i) == Atom::CARBON && nRingish <= 1) return 6;
  return 0;
}

/// BFS of the nucleobase from χ, never crossing back onto C1′.
static Iarray WalkBase(Graph const& g, int chi, int c1) {
  std::set<int> skip;
  if (c1 >= 0) skip.insert(c1);
  std::set<int> seen;
  seen.insert(chi);
  Iarray order(1, chi);
  Iarray stack(1, chi);
  size_t q = 0;
  while (q < stack.size()) {
    int i = stack[q++];
    Iarray kids;
    Iarray const& hn = g.HeavyNbr(i);
    for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
      if (seen.find(*j) != seen.end()) continue;
      if (skip.find(*j) != skip.end()) continue;
      kids.push_back(*j);
    }
    std::sort(kids.begin(), kids.end(), [&g](int a, int b) {
      int pa = BasePriority(g, a);
      int pb = BasePriority(g, b);
      if (pa != pb) return pa < pb;
      if (g.Name(a) != g.Name(b)) return g.Name(a) < g.Name(b);
      return a < b;
    });
    for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j) {
      seen.insert(*j);
      order.push_back(*j);
      stack.push_back(*j);
    }
  }
  return order;
}

/// Emit atom i and every still-unused neighbor subtree (used for 2′ substituents).
static void EmitSubtree(Graph const& g, int i, int parent,
                        std::vector<char>& used, Iarray& order)
{
  if (used[i]) return;
  used[i] = 1;
  order.push_back(i);
  Iarray kids;
  Iarray const& nbr = g.Nbr(i);
  for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
    if (*j != parent && !used[*j]) kids.push_back(*j);
  }
  std::sort(kids.begin(), kids.end(), [&g](int a, int b) {
    if (g.IsHydrogen(a) != g.IsHydrogen(b))
      return !g.IsHydrogen(a) && g.IsHydrogen(b);
    if (g.Elt(a) != g.Elt(b)) return g.Elt(a) < g.Elt(b);
    return g.Name(a) < g.Name(b);
  });
  for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
    EmitSubtree(g, *j, i, used, order);
}

/// Canonical nucleic-acid walk used when the user requests naorder.
/** P → OP (+H) → O5′ → C5′ → C4′ → O4′ → C1′ → base from χ →
  * C3′ → C2′ → 2′ substituents → O3′ → any leftover atoms.
  */
Iarray CanonicalWalk(Graph const& g, Scaffold const& sc) {
  std::vector<char> used(g.Natom(), 0);
  Iarray order;
  EmitIdx(sc.Get(Scaffold::P), used, order);
  for (Iarray::const_iterator op = sc.op_.begin(); op != sc.op_.end(); ++op)
    EmitHeavyAndH(g, *op, used, order);
  EmitHeavyAndH(g, sc.Get(Scaffold::O5p), used, order);
  EmitHeavyAndH(g, sc.Get(Scaffold::C5p), used, order);
  EmitHeavyAndH(g, sc.Get(Scaffold::C4p), used, order);
  EmitHeavyAndH(g, sc.Get(Scaffold::O4p), used, order);
  EmitHeavyAndH(g, sc.Get(Scaffold::C1p), used, order);
  int chi = sc.Get(Scaffold::CHI);
  if (chi >= 0) {
    Iarray base = WalkBase(g, chi, sc.Get(Scaffold::C1p));
    for (Iarray::const_iterator i = base.begin(); i != base.end(); ++i)
      EmitHeavyAndH(g, *i, used, order);
  }
  EmitHeavyAndH(g, sc.Get(Scaffold::C3p), used, order);
  int c2 = sc.Get(Scaffold::C2p);
  EmitHeavyAndH(g, c2, used, order);
  if (c2 >= 0) {
    std::set<int> skip;
    if (sc.Get(Scaffold::C1p) >= 0) skip.insert(sc.Get(Scaffold::C1p));
    if (sc.Get(Scaffold::C3p) >= 0) skip.insert(sc.Get(Scaffold::C3p));
    if (sc.Get(Scaffold::O4p) >= 0) skip.insert(sc.Get(Scaffold::O4p));
    Iarray kids;
    Iarray const& nbr = g.Nbr(c2);
    for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
      if (!used[*j] && skip.find(*j) == skip.end()) kids.push_back(*j);
    }
    std::sort(kids.begin(), kids.end(), [&g](int a, int b) {
      if (g.IsHydrogen(a) != g.IsHydrogen(b))
        return !g.IsHydrogen(a) && g.IsHydrogen(b);
      if (g.Elt(a) != g.Elt(b)) return g.Elt(a) < g.Elt(b);
      return g.Name(a) < g.Name(b);
    });
    for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
      EmitSubtree(g, *j, c2, used, order);
  }
  EmitHeavyAndH(g, sc.Get(Scaffold::O3p), used, order);
  for (int i = 0; i < g.Natom(); i++) {
    if (!used[i]) EmitIdx(i, used, order);
  }
  return order;
}

/// Count unused heavy atoms reachable from start without crossing parent or used.
static int HeavySubtreeSize(Graph const& g, int start, int parent,
                            std::vector<char> const& used)
{
  std::vector<char> seen(g.Natom(), 0);
  Iarray stack(1, start);
  seen[start] = 1;
  int count = 0;
  while (!stack.empty()) {
    int i = stack.back();
    stack.pop_back();
    if (used[i] && i != start) continue;
    if (!g.IsHydrogen(i)) count++;
    Iarray const& nbr = g.Nbr(i);
    for (Iarray::const_iterator j = nbr.begin(); j != nbr.end(); ++j) {
      if (*j == parent || seen[*j] || (used[*j] && *j != start)) continue;
      seen[*j] = 1;
      stack.push_back(*j);
    }
  }
  return count;
}

/// Walk a side-chain branch the way amino19.lib does: heavy, its hydrogens, then
/// child heavies (smaller unused subtree first, then name).
static void EmitAaBranch(Graph const& g, int i, int parent,
                         std::vector<char>& used, Iarray& order,
                         std::set<int> const& skip)
{
  if (i < 0 || used[i] || skip.find(i) != skip.end()) return;
  EmitHeavyAndH(g, i, used, order);
  Iarray kids;
  Iarray const& hn = g.HeavyNbr(i);
  for (Iarray::const_iterator j = hn.begin(); j != hn.end(); ++j) {
    if (*j == parent || used[*j] || skip.find(*j) != skip.end()) continue;
    kids.push_back(*j);
  }
  std::sort(kids.begin(), kids.end(), [&g, &used, i](int a, int b) {
    int sa = HeavySubtreeSize(g, a, i, used);
    int sb = HeavySubtreeSize(g, b, i, used);
    if (sa != sb) return sa < sb;
    return g.Name(a) < g.Name(b);
  });
  for (Iarray::const_iterator j = kids.begin(); j != kids.end(); ++j)
    EmitAaBranch(g, *j, i, used, order, skip);
}

/// Canonical amino-acid walk used when the user requests aaorder (ff19SB).
/** N (+ amide H) → [Pro: CD-ring to CB] → CA (+ HA) → CB side chain → C → O.
  * Leftovers (OXT, caps) append in index order.
  */
Iarray CanonicalAaWalk(Graph const& g, Scaffold const& sc) {
  std::vector<char> used(g.Natom(), 0);
  Iarray order;
  int n  = sc.Aa(Scaffold::AA_N);
  int ca = sc.Aa(Scaffold::AA_CA);
  int c  = sc.Aa(Scaffold::AA_C);
  int o  = sc.Aa(Scaffold::AA_O);
  int cb = sc.Aa(Scaffold::AA_CB);
  std::set<int> skip;
  if (n  >= 0) skip.insert(n);
  if (ca >= 0) skip.insert(ca);
  if (c  >= 0) skip.insert(c);
  if (o  >= 0) skip.insert(o);

  EmitHeavyAndH(g, n, used, order);
  if (sc.isPro_ && n >= 0) {
    Iarray const& nn = g.HeavyNbr(n);
    for (Iarray::const_iterator j = nn.begin(); j != nn.end(); ++j) {
      if (*j != ca && !used[*j])
        EmitAaBranch(g, *j, n, used, order, skip);
    }
  }
  EmitHeavyAndH(g, ca, used, order);
  if (!sc.isPro_) {
    if (cb >= 0)
      EmitAaBranch(g, cb, ca, used, order, skip);
    else if (ca >= 0) {
      Iarray const& can = g.HeavyNbr(ca);
      for (Iarray::const_iterator j = can.begin(); j != can.end(); ++j) {
        if (*j != n && *j != c && !used[*j] && !g.IsHydrogen(*j))
          EmitAaBranch(g, *j, ca, used, order, skip);
      }
    }
  }
  EmitHeavyAndH(g, c, used, order);
  EmitHeavyAndH(g, o, used, order);
  for (int i = 0; i < g.Natom(); i++) {
    if (!used[i]) EmitIdx(i, used, order);
  }
  return order;
}

/// Template order: file order, NA walk, or amino-acid walk (ff19SB).
Iarray ParentOrder(Graph const& tpl, Scaffold const& sc, bool useNaOrder, bool useAaOrder) {
  if (useNaOrder)
    return CanonicalWalk(tpl, sc);
  if (useAaOrder)
    return CanonicalAaWalk(tpl, sc);
  Iarray order(tpl.Natom());
  for (int i = 0; i < tpl.Natom(); i++)
    order[i] = i;
  return order;
}

} // namespace

/// Public wrapper: detect the scaffold of \p top and return CanonicalWalk.
int TemplateMatch::CanonicalNaOrder(Topology const& top, Iarray& order) const {
  Graph g(top);
  Scaffold sc;
  DetectScaffold(g, sc);
  order = CanonicalWalk(g, sc);
  return 0;
}

int TemplateMatch::CanonicalAaOrder(Topology const& top, Iarray& order) const {
  Graph g(top);
  Scaffold sc;
  DetectScaffold(g, sc);
  order = CanonicalAaWalk(g, sc);
  return 0;
}

/// Detect scaffolds, map graphs, build the output permutation, count partial-map stats.
int TemplateMatch::Match(Topology const& tgtTop, Topology const& tplTop, Result& out) const {
  out = Result();
  if (tgtTop.Natom() < 1 || tplTop.Natom() < 1) {
    mprinterr("Error: timap: empty topology.\n");
    return 1;
  }
  Graph tgt(tgtTop);
  Graph tpl(tplTop);
  Scaffold tgtSc, tplSc;
  DetectScaffold(tgt, tgtSc);
  DetectScaffold(tpl, tplSc);
  out.tgtKind_ = tgtSc.kind_;
  out.tplKind_ = tplSc.kind_;
  out.notes_.insert(out.notes_.end(), tgtSc.notes_.begin(), tgtSc.notes_.end());

  if (debug_ > 0) {
    mprintf("DEBUG: timap tgt kind=%s family=%s 2'=%s ok=%i\n",
            tgtSc.kind_.c_str(), tgtSc.family_.c_str(), tgtSc.twoPrime_.c_str(), (int)tgtSc.ok_);
    mprintf("DEBUG: timap tpl kind=%s family=%s 2'=%s ok=%i\n",
            tplSc.kind_.c_str(), tplSc.family_.c_str(), tplSc.twoPrime_.c_str(), (int)tplSc.ok_);
    for (int r = 0; r < Scaffold::NROLES; r++) {
      int ti = tgtSc.Get((Scaffold::Role)r);
      int pi = tplSc.Get((Scaffold::Role)r);
      mprintf("DEBUG:   %-4s  tgt=%s  tpl=%s\n", Scaffold::RoleStr((Scaffold::Role)r),
              ti < 0 ? "-" : tgt.Name(ti).c_str(),
              pi < 0 ? "-" : tpl.Name(pi).c_str());
    }
  }

  MapGraphs(tgt, tpl, tgtSc, tplSc, out.mapping_, seed_);
  Iarray parent = ParentOrder(tpl, tplSc, useNaOrder_, useAaOrder_);
  out.outputOrder_ = OutputOrder(tgt, tpl, tplSc, out.mapping_, parent, anchorName_);
  out.dual_ = DualOrder(tgt, tpl, tplSc, out.mapping_, parent, anchorName_);

  if ((int)out.outputOrder_.size() != tgt.Natom()) {
    mprinterr("Error: timap: output order size %zu != %i atoms.\n",
              out.outputOrder_.size(), tgt.Natom());
    return 1;
  }
  std::vector<char> seen(tgt.Natom(), 0);
  for (Iarray::const_iterator it = out.outputOrder_.begin(); it != out.outputOrder_.end(); ++it) {
    if (*it < 0 || *it >= tgt.Natom() || seen[*it]) {
      mprinterr("Error: timap: output order is not a permutation.\n");
      return 1;
    }
    seen[*it] = 1;
  }

  Iarray usedTpl(tpl.Natom(), 0);
  for (int n = 0; n < tgt.Natom(); n++) {
    if (out.mapping_[n] >= 0) {
      out.nMapped_++;
      usedTpl[out.mapping_[n]] = 1;
    } else {
      out.nInsertion_++;
    }
  }
  for (int r = 0; r < tpl.Natom(); r++) {
    if (!usedTpl[r]) out.nUnmappedTpl_++;
  }
  int expectDual = out.nMapped_ + out.nInsertion_ + out.nUnmappedTpl_;
  if ((int)out.dual_.size() != expectDual) {
    mprinterr("Error: timap: dual-topology size %zu != %i (mapped+ins+unmapped).\n",
              out.dual_.size(), expectDual);
    return 1;
  }
  return 0;
}
