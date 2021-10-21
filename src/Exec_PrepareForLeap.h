#ifndef INC_EXEC_PREPAREFORLEAP_H
#define INC_EXEC_PREPAREFORLEAP_H
#include "Exec.h"
#include "AtomMap.h"
#include <set>
#include <map>
#include <vector>
class AtomMask;
class CharMask;
class CpptrajFile;
/// Do common tasks to prepare a structure to be loaded into tleap 
class Exec_PrepareForLeap : public Exec {
  public:
    Exec_PrepareForLeap();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrepareForLeap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    /// Hold indices for sugar
    class Sugar;
    /// Hold indices for sugar link atoms
    class Link;

    typedef std::vector<int> Iarray;
    enum ResStatType { UNKNOWN = 0, VALIDATED, UNRECOGNIZED_SUGAR_LINKAGE, SUGAR_MISSING_C1X,
                       SUGAR_MISSING_RING_O };
    typedef std::vector<ResStatType> ResStatArray;

    inline void ChangeResName(Residue&, NameType const&) const;
    inline void ChangeAtomName(Atom&, NameType const&) const;

    /// Set a reduced PDB res to glycam map when dat file not found.
    void SetGlycamPdbResMap();
    /// Load PDB res to glycam map from dat file
    int LoadGlycamPdbResMap(std::string const&);
    /// Set PDB residue names recognized by Amber FFs
    void SetPdbResNames();
    /// Load PDB residue names recognized by Amber FFs from dat file
    int LoadPdbResNames(std::string const&);

    void LeapBond(int,int,Topology const&, CpptrajFile*) const;
//    int CalcStereocenterTorsion(double&, int, Topology const&, Frame const&) const;
    int FindRemainingChainCarbons(Iarray&, int, Topology const&, int, Iarray const&) const;
    /// Try to find any missing bonds to C1 atoms
//    int FindSugarC1Linkages(Sugar const&, Topology&, Frame const&) const;
    /// Determine orientation around anomeric carbon
    int CalcAnomericTorsion(double&, int, int, int, Iarray const&,
                            Topology const&, Frame const&) const;
    /// Determine orientation around anomeric reference carbon
    int CalcAnomericRefTorsion(double&, int, int, int, Iarray const&,
                               Topology const&, Frame const&) const;
    /// Determine orientation around configurational carbon
    int CalcConfigCarbonTorsion(double&, int, Iarray const&,
                                Topology const&, Frame const&) const;

    static int totalPriority(Topology const&, int, int, int, int, std::vector<bool>&);

    /// return type for the CalcChiralAtomTorsion routine
    enum ChiralRetType { ERR = 0, IS_S, IS_R };

    ChiralRetType CalcChiralAtomTorsion(double&, int, Topology const&, Frame const&) const;
    /// Error status for IdSugarRing
    enum IdSugarRingStatType { ID_OK = 0, ID_ERR, ID_MISSING_O };
    /// \return Sugar with atom indices set up
    Sugar IdSugarRing(int, Topology const&, IdSugarRingStatType&) const;
    int ChangePdbAtomNamesToGlycam(char, Residue const&, Topology&) const;

    /// Return type for DetermineAnomericForm
    enum AnomerRetType { A_ERR = 0, A_WARNING, IS_ALPHA, IS_BETA };
    /// \return Anomeric form of the sugar
    AnomerRetType DetermineAnomericForm(bool&, Sugar const&, Topology const&, Frame const&) const;
    /// \return Glycam linkage code for given link atoms
    std::string GlycamLinkageCode(std::set<Link> const&, Topology const&) const;
    /// Determine linkages for the sugar
    std::string DetermineSugarLinkages(Sugar const&, CharMask const&, Topology&, ResStatArray&,
                                       CpptrajFile*, std::set<BondType>&) const;
    /// Try to identify sugar name, form, and linkages
    int IdentifySugar(Sugar const&, Topology&, Frame const&, CharMask const&, CpptrajFile*, std::set<BondType>&);
    /// Try to find missing linkages to anomeric carbon in sugar.
    int FindSugarC1Linkages(int, int, Topology&, Frame const&) const;
    /// Determine if sugar has sulfates that need SO3 residue(s)
    int CheckForSugarSulfates(int, Iarray const&, Topology&, Frame&) const;
    /// Determine if sugar is terminal and need an ROH residue
    int CheckIfSugarIsTerminal(int, int, int, Topology&, Frame&) const;
    /// Attempt to fix any issues with sugars
    int FixSugarsStructure(std::string const&, Topology&, Frame&, bool, bool) const;

    int PrepareSugars(AtomMask&, Topology&, Frame const&, CpptrajFile*);
    int FindTerByBonds(Topology&, CharMask const&) const;
    int SearchForDisulfides(double, std::string const&, std::string const&, bool,
                            Topology&, Frame const&, CpptrajFile*);
    /// \return true if residue name is a recognized PDB name
    bool IsRecognizedPdbRes(NameType const&) const;
    /// \return Array of residue nums with unrecognized names
    Iarray GetUnrecognizedPdbResidues(Topology const&) const;
    /// \return Array indices of isolated unrecognized residues
    Iarray GetIsolatedUnrecognizedResidues(Topology const&, Iarray const&) const;
    /// Remove specified atoms
    int ModifyCoords(Topology&, Frame&, bool, char, std::string const&,
                     std::string const&, Iarray const&) const;
    int RemoveHydrogens(Topology&, Frame&) const;
    int DetermineHisProt(std::vector<NameType>&, Iarray&, Topology const&,
                         NameType const&, NameType const&,
                         NameType const&, NameType const&, NameType const&, NameType const&) const;


    typedef std::pair<NameType, char> PairType;
    typedef std::map<NameType, char> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to Glycam 1 char names

    typedef std::set<NameType> SetType;
    SetType pdb_res_names_; ///< PDB residue names recognized by Amber FFs

    ResStatArray resStat_;  ///< Contain status of each residue

    typedef std::pair<NameType, NameType> NamePairType;
    typedef std::map<NameType, NameType> NameMapType;
    /// Hold maps of pdb atom names to glycam atom names
    std::vector<NameMapType> pdb_glycam_name_maps_;

    typedef std::pair<char, int> ResIdxPairType;
    typedef std::map<char, int> ResIdxMapType;
    /// Map glycam residue chars to pdb-glycam atom name maps
    ResIdxMapType glycam_res_idx_map_;

    /// Map pdb residue names to glycam linkage residue names
    NameMapType pdb_glycam_linkageRes_map_;

    std::string leapunitname_;
    bool errorsAreFatal_;   ///< If false, try to skip errors.
    int debug_; ///< Debug level
    std::string solventResName_; ///< Solvent residue name
    std::string terminalHydroxylName_; ///< Terminal hydroxyl name
    AtomMap myMap_;
};
// ----- Sugar class ----------------------------------------------------------
class Exec_PrepareForLeap::Sugar {
  public:
    /// Base ring type
    enum RingTypeEnum { PYRANOSE = 0,  ///< Ring is 5 carbons, 1 oxygen
                        FURANOSE,      ///< Ring is 4 carbons, 1 oxygen
                        UNKNOWN_RING   ///< Some unknown ring type
                  };
    /// CONSTRUCTOR - residue number, incomplete setup
    Sugar(int);
    /// CONSTRUCTOR - res #, ring O, Anomeric, Anomeric Ref, Highest Sterocenter, ring atoms, chain atoms
    Sugar(int,int,int,int,int,Iarray const&,Iarray const&);

    int ResNum()               const { return rnum_; }
    int RingOxygenAtom()       const { return ring_oxygen_atom_; }
    int AnomericAtom()         const { return anomeric_atom_; }
    int AnomericRefAtom()      const { return ano_ref_atom_; }
    int HighestStereocenter()  const { return highest_stereocenter_; }
    RingTypeEnum RingType()    const { return ringType_; }
    Iarray const& RingAtoms()  const { return ring_atoms_; }
    int RingEndAtom()          const { return ring_atoms_.back(); }
    Iarray const& ChainAtoms() const { return chain_atoms_; }
//    typedef std::vector<int>::const_iterator const_iterator;
//    const_iterator ringbegin() const { return ring_atoms_.begin(); }
//    const_iterator ringend()   const { return ring_atoms_.end(); }

    bool NotSet() const { return (ring_oxygen_atom_ == -1); }
    /// \return Number of ring atoms
    unsigned int NumRingAtoms() const;
    void PrintInfo(Topology const&) const;
  private:
    int rnum_;                 ///< Residue index
    int ring_oxygen_atom_;     ///< Index of the ring oxygen atom
    int anomeric_atom_;        ///< Index of the anomeric C atom (ring start)
    int ano_ref_atom_;         ///< Index of the anomeric reference C atom
    int highest_stereocenter_; ///< Index of the highest stereocenter in the carbon chain
    RingTypeEnum ringType_;    ///< Will be set to ring type
    Iarray ring_atoms_;        ///< Index of all non-oxygen ring atoms
    Iarray chain_atoms_;       ///< Index of all chain carbon atoms (from anomeric carbon).
};
// ----- Link class ------------------------------------------------------------
class Exec_PrepareForLeap::Link {
  public:
    /// CONSTRUCTOR
    Link() : idx_(-1), position_(-1) {}
    /// CONSTRUCTOR - index, position (starting from 1 at the anomeric carbon)
    Link(int i, int p) : idx_(i), position_(p) {}
    /// \return Index in topology
    int Idx() const { return idx_; }
    /// \return Index in carbon chain (starting from 1 at the anomeric carbon)
    int Position() const { return position_; }
    /// First sort by position, then absolute index
    bool operator<(Link const& rhs) const {
      if (position_ == rhs.position_) {
        return (idx_ < rhs.idx_);
      } else {
        return position_ < rhs.position_;
      }
    }
  private:
    int idx_;      ///< Atom index in topology
    int position_; ///< Position in sugar carbon chain, starting from 1 at the anomeric carbon
};
#endif
