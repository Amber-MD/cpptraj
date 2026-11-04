#ifndef INC_EXEC_FIXATOMNAMES_H
#define INC_EXEC_FIXATOMNAMES_H
#include "Exec.h"
//#incl ude <map>
/// Fix topology so all atom names in residues are unique.
class Exec_FixAtomNames : public Exec {
  public:
    Exec_FixAtomNames() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_FixAtomNames(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
/*    typedef std::pair<NameType, NameType> NamePair;
    typedef std::map<NameType, NameType> NameMap;
    typedef std::pair<NameType, NameMap> ResPair;
    typedef std::map<NameType, NameMap> ResMap;

    ResMap resMap_; ///< Hold a map of residues that have had their atom names changed*/
    typedef std::vector<int> Iarray;
    typedef std::pair<NameType, Iarray> DuplicatePair;
    typedef std::map<NameType, Iarray> DuplicateMap;

    /// Hold a part of an atom name
    class NameToken {
      public:
        enum Type { STRING = 0, ///< A letter (A-Z, a-z)
                    INT,        ///< A digit (0-9)
                    OTHER,      ///< Any other character. Will not change.
                    UNKNOWN };
        NameToken(std::string const& s) : ival_(-1), sval_(s), type_(STRING) {}
        NameToken(int i) : ival_(i), sval_(""), type_(INT) {}
        NameToken(Type t, std::string const& s) : ival_(-1), sval_(s), type_(t) {}

        int Ival() const { return ival_; }
        std::string const& Sval() const { return sval_; }
        Type TokenType() const { return type_; }

        void SetIval(int i) { ival_ = i; }
      private:
        int ival_;
        std::string sval_;
        Type type_;
    };
    /// Hold all parts that make up an atom name
    class AtomName {
        typedef std::vector<NameToken> NameTokenArray;
      public:
        AtomName() {}
        //void AddToken(NameToken const& n) { tokens_.push_back( n ); }
        void PrintName() const;
        void PrintName(const char*) const;
        std::string NameString() const;
        int DissectName(NameType const&);
        int NameWidth() const;
        /// Naively increment the name
        int IncrementName(int);
        /// \return Terminal number, or -1 if no terminal number
        int TerminalNumber() const {
          NameTokenArray::const_reverse_iterator rit1 = tokens_.rbegin();
          while (rit1 != tokens_.rend()) {
            if (rit1->TokenType() == NameToken::INT) break;
            ++rit1;
          }
          if (rit1 == tokens_.rend()) return -1;
          return rit1->Ival();
        }
        /// \return true if name contains a prime (')
        bool HasPrime() const {
          for (NameTokenArray::const_iterator it1 = tokens_.begin(); it1 != tokens_.end(); ++it1)
            if (it1->TokenType() == NameToken::OTHER && it1->Sval() == "'") return true;
          return false;
        }
        /// Remove any leading digits
        int RemoveLeadingDigit();
        /// Set or add a terminal number
        int SetTerminalNumber(int, int);
        /// \return True if both atoms have a terminal number and the numbers match
        bool TerminalNumbersMatch(AtomName const& rhs) const {
          NameTokenArray::const_reverse_iterator rit1 = tokens_.rbegin();
          while (rit1 != tokens_.rend()) {
            if (rit1->TokenType() == NameToken::INT) break;
            ++rit1;
          }
          if (rit1 == tokens_.rend()) return false;
          NameTokenArray::const_reverse_iterator rit2 = rhs.tokens_.rbegin();
          while (rit2 != rhs.tokens_.rend()) {
            if (rit2->TokenType() == NameToken::INT) break;
            ++rit2;
          }
          if (rit2 == rhs.tokens_.rend()) return false;
          // At this point we know both names have a terminal digit.
          return (rit1->Ival() == rit2->Ival());
        }
        /// If both atoms have a terminal number, make this number match the given number
        bool MakeTerminalNumbersMatch(AtomName const& rhs) {
          NameTokenArray::reverse_iterator rit1 = tokens_.rbegin();
          while (rit1 != tokens_.rend()) {
            if (rit1->TokenType() == NameToken::INT) break;
            ++rit1;
          }
          if (rit1 == tokens_.rend()) return false;
          NameTokenArray::const_reverse_iterator rit2 = rhs.tokens_.rbegin();
          while (rit2 != rhs.tokens_.rend()) {
            if (rit2->TokenType() == NameToken::INT) break;
            ++rit2;
          }
          if (rit2 == rhs.tokens_.rend()) return false;
          // At this point we know both names have a terminal digit.
          rit1->SetIval( rit2->Ival() );
          return true;
        }
      private:
        NameTokenArray tokens_;
    };
    /// Hold info for a duplicated atom
    class DupAtom {
      public:
        DupAtom(int i, int s) : idx_(i), sum_(s) {}
        int Idx() const { return idx_; }
        int Sum() const { return sum_; }
        bool operator<(const DupAtom& rhs) const {
          if (sum_ == rhs.sum_)
            return (idx_ < rhs.idx_);
          else
            return (sum_ > rhs.sum_);
        }
      private:
        int idx_; ///< Atom index in topology
        int sum_; ///< Used to prioritize duplicated atoms; sum of atomic #s of bonded atoms
    };
    /// Used for sorting the duplicated name map
    class DupName {
      public:
        DupName(NameType const& n, Iarray const& i, int w) : weight_(w), aname_(n), idxs_(i) {}
        int weight_; ///< 9999 for C, atomic number for others
        NameType aname_;
        Iarray idxs_; ///< Atom indices in topology

        bool operator<(const DupName& rhs) const {
          if (weight_ == rhs.weight_) {
            if (aname_ == rhs.aname_) {
              return (idxs_.front() < rhs.idxs_.front());
            } else {
              return (aname_ < rhs.aname_);
            }
          } else {
            return (weight_ > rhs.weight_);
          }
        }
    };
        
    int change_aname(Topology&, int, NameType const&);

    static inline NameToken::Type get_type(char);

    int findDuplicatesInResidue(Topology const&, int);

    int fixDuplicateByIncrement(Topology&, int);

    int fixDuplicateSet(Topology&, Iarray const&, int);

    int fixChainedAtoms(Topology&, int, int);

    int RenameBondedHydrogens(Topology&, int);

    int MakeTerminalNumberMatch(Topology&, int, int, AtomName&, AtomName const&);

    int RenameAllHydrogens(Topology&, int);

    static const int maxNameWidth_;
    DuplicateMap duplicateMap_; ///< Hold atom names and indices for the current residue being fixed
};
#endif
