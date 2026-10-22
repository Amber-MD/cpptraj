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
      private:
        int ival_;
        std::string sval_;
        Type type_;
    };
    /// Hold all parts that make up an atom name
    class AtomName {
      public:
        AtomName() {}
        //void AddToken(NameToken const& n) { tokens_.push_back( n ); }
        void PrintName() const;
        int DissectName(NameType const&);
        int NameWidth() const;
        /// Naively increment the name
        int IncrementName(int);
      private:
        typedef std::vector<NameToken> NameTokenArray;
        NameTokenArray tokens_;
    };

    static inline NameToken::Type get_type(char);

    int fixDuplicate(Topology&, int, int) const;

};
#endif
