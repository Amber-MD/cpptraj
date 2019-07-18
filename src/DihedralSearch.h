#ifndef INC_DIHEDRALSEARCH_H
#define INC_DIHEDRALSEARCH_H
#include "MetaData.h"
#include "AtomMask.h"
// Forward declarations
class ArgList;
class Range;
class Topology;
/// Class that can be used to search for dihedral angles in a range.
// Thanks to C. Bergonzo for the NA angle definitions.
class DihedralSearch {
    class DihedralMask;
  public:
    /// Const iterator over found dihedrals (for frame processing)
    typedef std::vector<DihedralMask>::const_iterator mask_it;
    mask_it begin()  const { return dihedrals_.begin();     }
    mask_it end()    const { return dihedrals_.end();       }
    /// Number of found dihedrals
    int Ndihedrals() const { return (int)dihedrals_.size(); }
    /// Recognized dihedral types
    typedef MetaData::scalarType DihedralType;
    DihedralSearch();
    /// COPY CONSTRUCTOR - Set up for same types as input
    DihedralSearch(DihedralSearch const&);
    DihedralSearch& operator=(DihedralSearch const&);
    static void ListKnownTypes();
    static void OffsetHelp();
    static DihedralType GetType(std::string const&);
    /// Add a known dihedral type to be searched for.
    int SearchFor(DihedralType);
    /// Search for known dihedral type keywords.
    void SearchForArgs(ArgList&);
    /// Contains keywords for SearchNewTypeArgs()
    static const char* newTypeArgsHelp_;
    /// Search for new type via args
    int SearchForNewTypeArgs(ArgList&);
    /// Add a new dihedral type to be searched for.
    int SearchForNewType(int, std::string const&, std::string const&, std::string const&, 
                         std::string const&, std::string const&);
    /// Add all dihedral types if none have been added yet.
    int SearchForAll();
    /// \return True if no dihedral tokens
    bool NoDihedralTokens() const { return dihedralTokens_.empty(); }
    /// Find specified dihedrals for residues in Range.
    int FindDihedrals(Topology const&, Range const&);
    /// Clear found dihedrals and tokens.
    void Clear();
    /// Clear found dihedrals only.
    void ClearFound() { dihedrals_.clear(); }
    /// Print dihedrals currently being searched for.
    void PrintTypes();
    /// \return Mask of atoms that will move upon rotation.
    static AtomMask MovingAtoms(Topology const&, int, int);
  private:
    static const DihedralType D_FIRST;
    static const DihedralType D_END;
    struct DIH_TYPE;
    static const DIH_TYPE DIH[];
    class DihedralToken;
    std::vector<DihedralToken> dihedralTokens_; ///< Dihedrals to search for
    std::vector<DihedralMask> dihedrals_;       ///< Contains atom #s for each found dihedral
};
// ----- PRIVATE CLASS HEADERS -------------------------------------------------
/// Hold dihedral type information used for searching.
class DihedralSearch::DihedralToken {
  public:
    DihedralToken() : centerIdx_(2), type_(MetaData::UNDEFINED) {}
    /// Constructor for custom dihedral type
    DihedralToken(int, NameType const&, NameType const&, NameType const&, NameType const&,
                  std::string const&);
    /// Constructor for default dihedral type
    DihedralToken(DIH_TYPE const&);
    /// \return DihedralMask with 4 atoms corresponding to this dihedral for given residue.
    DihedralMask FindDihedralAtoms(Topology const&, int) const;
    std::string const& Name() const { return name_; }
    DihedralType Type()       const { return type_; }
    void SetAtomName(int i, NameType const& n) { aname_[i] = n; }
  private:
    int centerIdx_;           ///< Index of the "central" dihedral atom (determines res#).
    NameType aname_[4];       ///< Dihedral atom names. 
    std::string name_;        ///< Dihedral type name.
    DihedralType type_;       ///< Dihedral type.
};
/// Hold dihedral atom #s, residue #, and type name.
class DihedralSearch::DihedralMask {
  public:
    DihedralMask();
    DihedralMask(int,int,int,int,int,std::string const&,DihedralType);
    int A0()                  const { return a0_;         }
    int A1()                  const { return a1_;         }
    int A2()                  const { return a2_;         }
    int A3()                  const { return a3_;         }
    int ResNum()              const { return res_;        }
    std::string const& Name() const { return name_;       }
    bool None()               const { return (a0_ == -1); }
    DihedralType Type()       const { return type_;       }
  private:
    int a0_, a1_, a2_, a3_, res_;
    std::string name_;
    DihedralType type_; ///< For mapping to types in DataSet, statistics analysis
};
#endif
