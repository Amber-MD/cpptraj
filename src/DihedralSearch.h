#ifndef INC_DIHEDRALSEARCH_H
#define INC_DIHEDRALSEARCH_H
#include "Topology.h"
#include "Range.h"
#include "ArgList.h"
/// Class that can be used to search for dihedral angles in a range.
// Thanks to C. Bergonzo for the NA angle definitions.
class DihedralSearch {
    class DihedralMask;
  public:
    typedef std::vector<DihedralMask>::const_iterator mask_it;
    mask_it begin()  const { return dihedrals_.begin();     }
    mask_it end()    const { return dihedrals_.end();       }
    int Ndihedrals() const { return (int)dihedrals_.size(); }
    /// Recognized dihedral types
    enum DihedralType { PHI = 0, PSI, CHIP, OMEGA, ALPHA, BETA, GAMMA, DELTA,
                        EPSILON, ZETA, NU1, NU2, CHIN, NDIHTYPE };
    DihedralSearch();
    static void ListKnownTypes();
    static void OffsetHelp();
    static DihedralType GetType(std::string const&);
    /// Add a known dihedral type to be searched for.
    int SearchFor(DihedralType);
    /// Search for known dihedral type keywords.
    void SearchForArgs(ArgList&);
    /// Add a new dihedral type to be searched for.
    int SearchForNewType(int, std::string const&, std::string const&, std::string const&, 
                         std::string const&, std::string const&);
    /// Add all dihedral types if none have been added yet.
    int SearchForAll();
    /// Find specified dihedrals for residues in Range.
    int FindDihedrals(Topology const&, Range const&);
    /// Clear found dihedrals and tokens.
    void Clear();
    /// Clear found dihedrals only.
    void ClearFound();
    /// Print dihedrals currently being searched for.
    void PrintTypes();
    /// \return Mask of atoms that will move upon rotation.
    static AtomMask MovingAtoms(Topology const&, int, int);
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
        DihedralToken(int, NameType const&, NameType const&, NameType const&, NameType const&,
                      std::string const&);
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
