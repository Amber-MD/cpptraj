#ifndef INC_TIMATCH_H
#define INC_TIMATCH_H
#include <string>
#include <vector>
class Topology;
/// Map a target topology onto a user-supplied template (engine for the timap command).
/** Source encoding is UTF-8. Unicode (O5′, χ, λ, 2′) appears in comments only;
  * string literals in the implementation stay ASCII so the file compiles
  * without a compiler UTF-8 flag.
  *
  * The template atom order *is* the thermodynamic-integration (TI) shared-atom
  * order. Target atoms that correspond to a template atom occupy that same
  * slot. Extra atoms on the target — a 2′-OH on RNA vs DNA, a phenolic OH,
  * selenium in place of sulfur — are insertions. Insertions that are bonded
  * to a mapped parent are emitted immediately after that parent. Completely
  * unmatched leftovers are placed just before an anchor atom (default O3′).
  *
  * Matching is graph-based and is expected to be partial (unlike atommap,
  * which requires a 1:1 correspondence):
  *   1. Seed by unique atom names, nucleic-acid scaffold roles, or amino-acid
  *      N / CA / C / O / CB (ff19SB amino19.lib baseline).
  *   2. Grow by pairing unmatched neighbors that share a unique signature
  *      (element, heavy degree, sorted neighbor elements).
  *   3. Attach hydrogens to already-mapped heavy atoms, ordered by name.
  *
  * The same engine is used for nucleotides, amino acids, and small molecules.
  * Nucleic-acid scaffold detection is opportunistic: if no furanose / χ is
  * found, the residue is checked for a peptide N–CA–C=O motif (kind "amino");
  * otherwise kind is "unknown" and unique-name seeding is used.
  *
  * Implementation details, invariants, and "where to change what" live in the
  * file-level comment at the top of TIMatch.cpp - start there when maintaining.
  */
class TIMatch {
  public:
    typedef std::vector<int> Iarray;
    /// How the correspondence is started before the grow pass.
    enum SeedType { SEED_AUTO = 0, SEED_NAMES, SEED_NA, SEED_AA, SEED_NONE };

    /// Outcome of Match(): correspondence plus the permutation of the target.
    class Result {
      public:
        Result() : nMapped_(0), nInsertion_(0), nUnmappedTpl_(0) {}
        Iarray mapping_;     ///< tgt index → template index; −1 = insertion
        Iarray outputOrder_; ///< Map[newatom] = old tgt atom (full permutation)
        /// Dual-topology slots: tpl-only (dummy on λ=1), tgt-only (dummy on λ=0), or shared.
        struct DualSlot {
          DualSlot() : tpl_(-1), tgt_(-1) {}
          DualSlot(int r, int n) : tpl_(r), tgt_(n) {}
          int tpl_; ///< template atom, or −1 if this slot is a target insertion
          int tgt_; ///< target atom, or −1 if this slot is an unmatched template atom
          bool IsShared()  const { return tpl_ >= 0 && tgt_ >= 0; }
          bool IsTplOnly() const { return tpl_ >= 0 && tgt_ <  0; }
          bool IsTgtOnly() const { return tpl_ <  0 && tgt_ >= 0; }
        };
        std::vector<DualSlot> dual_;
        int nMapped_;        ///< target atoms that correspond to a template atom
        int nInsertion_;     ///< target atoms with no template partner
        int nUnmappedTpl_;   ///< template atoms with no target partner
        std::string tgtKind_; ///< nucleotide | sugar | base | amino | unknown
        std::string tplKind_;
        std::vector<std::string> notes_; ///< scaffold warnings (empty sugar, …)
    };

    TIMatch();

    void SetSeed(SeedType s)                 { seed_ = s; }
    void SetDebug(int d)                     { debug_ = d; }
    void SetUseNaOrder(bool b)               { useNaOrder_ = b; }
    void SetUseAaOrder(bool b)               { useAaOrder_ = b; }
    void SetAnchorName(std::string const& n) { anchorName_ = n; }

    /// Map tgt onto tpl. Partial maps are success (return 0) as long as the
    /// output order is a permutation of the target.
    int Match(Topology const& tgt, Topology const& tpl, Result& out) const;
    /// Nucleic-acid canonical walk of top: P → OP → O5′ → C5′ → C4′ → O4′ →
    /// C1′ → base from χ → C3′ → C2′ (+ 2′ substituents) → O3′.
    int CanonicalNaOrder(Topology const& top, Iarray& order) const;
    /// Amino-acid canonical walk of top (ff19SB amino19.lib): N → H → CA → HA →
    /// side chain from CB (each heavy, then its hydrogens) → C → O.
    /// Proline walks N → CD → … → CB → CA → C → O.
    int CanonicalAaOrder(Topology const& top, Iarray& order) const;

    static const char* SeedStr(SeedType);
    static SeedType SeedFromString(std::string const&);

  private:
    SeedType seed_;
    int debug_;
    bool useNaOrder_;
    bool useAaOrder_;
    std::string anchorName_;
};
#endif

