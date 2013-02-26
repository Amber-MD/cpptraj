#ifndef INC_DIHEDRALSEARCH_H
#define INC_DIHEDRALSEARCH_H
#include "Topology.h"
#include "Range.h"
class DihedralSearch {
    class DihedralMask;
  public:
    typedef std::vector<DihedralMask>::const_iterator mask_it;
    mask_it begin()  const { return dihedrals_.begin();     }
    mask_it end()    const { return dihedrals_.end();       }
    int Ndihedrals() const { return (int)dihedrals_.size(); }
    /// Recognized dihedral types
    enum DihedralType { PHI = 0, PSI, NDIHTYPE };
    DihedralSearch();
    static void ListKnownTypes();
    static DihedralType GetType(std::string const&);
    /// Add a new dihedral type to be searched for.
    int SearchFor(DihedralType);
    /// Add all dihedral types if none have been added yet.
    int SearchForAll();
    /// Find specified dihedrals for residues in Range.
    int FindDihedrals(Topology const&, Range const&);
    /// Clear dihedrals
    void Clear();
    /// Print dihedrals currently being searched for.
    void PrintTypes();
  private:
    /// Hold dihedral atom #s, residue #, and type name.
    class DihedralMask {
      public:
        DihedralMask();
        DihedralMask(int,int,int,int,int,std::string const&);
        int A0()                  const { return a0_;         }
        int A1()                  const { return a1_;         }
        int A2()                  const { return a2_;         }
        int A3()                  const { return a3_;         }
        int ResNum()              const { return res_;        }
        std::string const& Name() const { return name_;       }
        bool None()               const { return (a0_ == -1); }
      private:
        int a0_, a1_, a2_, a3_, res_;
        std::string name_;
    };
    /// Hold dihedral type information used for searching.
    class DihedralToken {
      public:
        DihedralToken();
        DihedralToken(int, const char*, const char*, const char*, const char*, const char*);
        /// \return mask with 4 atoms corresponding to dihedral for specified residue.
        DihedralMask FindDihedralAtoms(Topology const&, int);
        std::string const& Name() { return name_; }
      private:
        int offset_;       ///< -1|0|1: Dihedral starts at prev.|stays in current|ends at next res.
        NameType aname0_;  ///< Dihedral 1st atom name.
        NameType aname1_;  ///< Dihedral 2nd atom name.
        NameType aname2_;  ///< Dihedral 3rd atom name.
        NameType aname3_;  ///< Dihedral 4th atom name.
        std::string name_; ///< Dihedral name
    };
    std::vector<DihedralToken> dihedralTokens_; ///< Dihedrals to search for
    std::vector<DihedralMask> dihedrals_;       ///< Contains atom #s for each found dihedral
};
#endif
