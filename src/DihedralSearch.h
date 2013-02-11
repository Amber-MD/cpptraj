#ifndef INC_DIHEDRALSEARCH_H
#define INC_DIHEDRALSEARCH_H
#include "Topology.h"
#include "Range.h"
class DihedralSearch {
  public:
    typedef std::vector<AtomMask>::const_iterator mask_it;
    mask_it begin() { return dihedrals_.begin(); }
    mask_it end()   { return dihedrals_.end();   }
    typedef std::vector< std::pair<int,std::string> >::const_iterator dihres_it;
    dihres_it dihbegin() { return dihRes_.begin(); }
    dihres_it dihend()   { return dihRes_.end();   }
    /// Recognized dihedral types
    enum DihedralType { PHI = 0, PSI, NDIHTYPE };
    DihedralSearch();
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
    class DihedralToken {
      public:
        DihedralToken();
        DihedralToken(int, const char*, const char*, const char*, const char*, const char*);
        /// \return mask with 4 atoms corresponding to dihedral for specified residue.
        AtomMask FindDihedralAtoms(Topology const&, int);
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
    std::vector<AtomMask> dihedrals_;           ///< Contains atom #s for each dihedral
    ///< Contain residue # and dihedral name for each dihedral.
    std::vector< std::pair<int,std::string> > dihRes_;
};
#endif
