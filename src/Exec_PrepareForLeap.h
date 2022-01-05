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
    /// Hold information for the various functional group types (FunctionalGroupType)
    class FunctionalGroup;
    /// Hold general info for a specific sugar
    class SugarToken;

    typedef std::vector<int> Iarray;
    enum FunctionalGroupType { G_SO3 = 0, G_CH3, G_ACX, G_OH, G_OME, UNRECOGNIZED_GROUP };
    enum ResStatType { UNKNOWN = 0,
                       VALIDATED,
                       SUGAR_UNRECOGNIZED_LINK_RES,
                       SUGAR_UNRECOGNIZED_LINKAGE,
                       SUGAR_NO_LINKAGE,
                       SUGAR_NO_CHAIN_FOR_LINK,
                       SUGAR_NAME_MISMATCH,
                       //SUGAR_MISSING_C1X,
                       SUGAR_SETUP_FAILED };
    typedef std::vector<ResStatType> ResStatArray;

    /// Base ring type
    enum RingTypeEnum { PYRANOSE = 0,  ///< Ring is 5 carbons, 1 oxygen
                        FURANOSE,      ///< Ring is 4 carbons, 1 oxygen
                        UNKNOWN_RING   ///< Some unknown ring type
                  };
    /// Sugar form
    enum FormTypeEnum { ALPHA = 0, BETA, UNKNOWN_FORM };
    /// Sugar chirality
    enum ChirTypeEnum { IS_D, IS_L, UNKNOWN_CHIR };

    /// Keep synced with FunctionalGroupType
    static const char* FunctionalGroupStr_[];
    /// Keep synced with RingTypeEnum
    static const char* ringstr_[];
    /// Keep synced with FormTypeEnum
    static const char* formstr_[];
    /// Keep synced with ChirTypeEnum
    static const char* chirstr_[];

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

    /// \return Sugar with atom indices set up
    Sugar IdSugarRing(int, Topology const&) const;
    /// Change PDB atom names to Glycam names
    int ChangePdbAtomNamesToGlycam(char, Residue const&, Topology&, FormTypeEnum) const;
    /// Determine form/chirality for furanose
    int DetermineUpOrDown(SugarToken&, Sugar const&, Topology const&, Frame const&) const;
    /// Determine form/chirliaty for pyranose 
    int DetermineAnomericForm(SugarToken&, Sugar&, Topology const&, Frame const&) const;
    /// \return Glycam linkage code for given link atoms
    std::string GlycamLinkageCode(std::set<Link> const&, Topology const&) const;
    /// Determine linkages for the sugar
    std::string DetermineSugarLinkages(Sugar const&, CharMask const&, Topology&, ResStatArray&,
                                       CpptrajFile*, std::set<BondType>&) const;
    /// Try to identify sugar name, form, and linkages
    int IdentifySugar(Sugar&, Topology&, Frame const&, CharMask const&, CpptrajFile*, std::set<BondType>&);
    /// Try to find missing linkages to anomeric carbon in sugar.
    int FindSugarC1Linkages(int, int, Topology&, Frame const&) const;
    /// \return identity of the group bonded to given atom
    FunctionalGroupType IdFunctionalGroup_Silent(Iarray&, int, int, int, Topology const&) const;
    /// \return identity of the group bonded to given atom, print to stdout
    FunctionalGroupType IdFunctionalGroup(Iarray&, int, int, int, Topology const&) const;
    /// Determine if sugar has sulfates that need SO3 residue(s)
    int CheckForFunctionalGroups(Sugar&, Topology&, Frame&) const;
    /// Determine if sugar is terminal and need an ROH residue
    int CheckIfSugarIsTerminal(Sugar&, Topology&, Frame&) const;
    /// Attempt to fix any issues with sugars
    int FixSugarsStructure(std::vector<Sugar>&, std::string const&, Topology&, Frame&,
                           bool, bool) const;

    int PrepareSugars(std::string const&, std::vector<Sugar>&, Topology&, Frame const&, CpptrajFile*);
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
    /** Try to determine protonation state of histidines from any hydrogens present. */
    int DetermineHisProt(Topology&,
                         NameType const&, NameType const&,
                         NameType const&, NameType const&, NameType const&, NameType const&) const;
    /// Run leap to generate topology, perform any modifications
    int RunLeap(std::string const&, std::string const&) const;

    typedef std::pair<NameType, SugarToken> PairType;
    typedef std::map<NameType, SugarToken> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to sugar information tokens

    typedef std::set<NameType> SetType;
    SetType pdb_res_names_; ///< PDB residue names recognized by Amber FFs

    ResStatArray resStat_;  ///< Contain status of each residue

    typedef std::pair<NameType, NameType> NamePairType;
    typedef std::map<NameType, NameType> NameMapType;
  
    static void PrintAtomNameMap(const char*, std::vector<NameMapType> const&);

    /// Hold maps of pdb atom names to glycam atom names; multiple residues may share a map
    std::vector<NameMapType> pdb_glycam_name_maps_;
    /// Hold maps of pdb atom names to glycam atom names (res in alpha form)
    std::vector<NameMapType> pdb_glycam_name_maps_A_;
    /// Hold maps of pdb atom names to glycam atom names (res in beta form)
    std::vector<NameMapType> pdb_glycam_name_maps_B_;

    typedef std::pair<char, int> ResIdxPairType;
    typedef std::map<char, int> ResIdxMapType;
    /// Map glycam residue chars to pdb-glycam atom name map index (pdb_glycam_name_maps_)
    ResIdxMapType glycam_res_idx_map_;

    /// Map pdb residue names to glycam linkage residue names
    NameMapType pdb_glycam_linkageRes_map_;

    std::string leapunitname_;
    bool errorsAreFatal_;   ///< If false, try to skip errors.
    bool hasGlycam_;        ///< If true, assume sugars already have glycam names
    int debug_;             ///< Debug level
    std::string solventResName_; ///< Solvent residue name
    std::string terminalHydroxylName_; ///< Terminal hydroxyl name
    AtomMap myMap_;
    std::vector<FunctionalGroup> functionalGroups_; ///< Recognized functional groups (FunctionalGroupType).
};
// ----- Sugar class ----------------------------------------------------------
class Exec_PrepareForLeap::Sugar {
  public:

    /// Sugar status. Keep synced with StatTypeStr_
    enum StatType { SETUP_OK = 0,    ///< Regular sugar, setup complete
                    MISSING_O,       ///< Could not find ring oxygen
                    MULTIPLE_O,      ///< Multiple potential ring oxygen atoms
                    MISSING_CHAIN,   ///< Could not find all chain carbons
                    MISSING_ANO_REF, ///< Missing anomeric reference atom
                    MISSING_CONFIG,  ///< Missing configurational carbon
                    MISSING_C1X      ///< Missing C1 X substituent
                  };
    /// CONSTRUCTOR - Status and residue first atom, incomplete setup
    Sugar(StatType, int);
    /// CONSTRUCTOR - Status, ring O, anomeric, ring atoms, chain atoms, incomplete setup
    Sugar(StatType, int, int, Iarray const&,Iarray const&);
    /// CONSTRUCTOR - ring O, Anomeric, Anomeric Ref, Highest Sterocenter, ring atoms, chain atoms, isMissingAtoms
    Sugar(int,int,int,int,Iarray const&,Iarray const&, bool);

    inline int ResNum(Topology const&) const;
    StatType Status()          const { return stat_; }
    int RingOxygenAtom()       const { return ring_oxygen_atom_; }
    int AnomericAtom()         const { return anomeric_atom_; }
    int AnomericRefAtom()      const { return ano_ref_atom_; }
    int HighestStereocenter()  const { return highest_stereocenter_; }
    RingTypeEnum RingType()    const { return ringType_; }
    bool IsMissingAtoms()      const { return isMissingAtoms_; }
    Iarray const& RingAtoms()  const { return ring_atoms_; }
    int RingEndAtom()          const { return ring_atoms_.back(); }
    Iarray const& ChainAtoms() const { return chain_atoms_; }
//    typedef std::vector<int>::const_iterator const_iterator;
//    const_iterator ringbegin() const { return ring_atoms_.begin(); }
//    const_iterator ringend()   const { return ring_atoms_.end(); }

    bool NotSet() const { return (stat_ != SETUP_OK); }
    /// \return Number of ring atoms
    unsigned int NumRingAtoms() const;
    void PrintInfo(Topology const&) const;
    /// Remap internal indices according to given atom map.
    void RemapIndices(Iarray const&, int, int);

    void SetStatus(StatType s) { stat_ = s; }
  private:
    /// Strings corresponding to StatType
    static const char* StatTypeStr_[];

    StatType stat_;            ///< Setup status
    int ring_oxygen_atom_;     ///< Index of the ring oxygen atom
    int anomeric_atom_;        ///< Index of the anomeric C atom (ring start)
    int ano_ref_atom_;         ///< Index of the anomeric reference C atom
    int highest_stereocenter_; ///< Index of the highest stereocenter in the carbon chain
    RingTypeEnum ringType_;    ///< Will be set to ring type
    bool isMissingAtoms_;      ///< True if original PDB indicates sugar has missing atoms.
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
// ----- FunctionalGroup class -------------------------------------------------
class Exec_PrepareForLeap::FunctionalGroup {
  public:
    FunctionalGroup();
  private:
    NameType resname_;                   ///< Functional group residue name.
    std::vector<NameType> anames_;       ///< Functional group atom names. Heavy atoms first.
    Atom::AtomicElementType chargeAtom_; ///< Element of atom which needs charge adjusted.
    double chargeOffset_;                ///< Charge offset for adjusting charge.
};
// ----- SugarToken class ------------------------------------------------------
class Exec_PrepareForLeap::SugarToken {
  public:
    /// CONSTRUCTOR
    SugarToken();
    /// CONSTRUCTOR - name, glycam code, form, chirality, ring type
    SugarToken(std::string const&, std::string const&, FormTypeEnum, ChirTypeEnum, RingTypeEnum);
    /// CONSTRUCTOR - ring type
    SugarToken(RingTypeEnum);
    /// /return <res>, set up from line: '<res> <code> <form> <chir> <ring> <name>'
    std::string SetFromLine(ArgList const&);
    /// Print token info to stdout
    void PrintInfo(std::string const&) const;

    std::string const& FullName()   const { return name_; }
    std::string const& GlycamCode() const { return glycamCode_; }
    FormTypeEnum Form()             const { return form_; }
    ChirTypeEnum Chirality()        const { return chir_; }
    RingTypeEnum RingType()         const { return ring_; }

    void SetChirality(ChirTypeEnum c) { chir_ = c; }
    void SetForm(FormTypeEnum f)      { form_ = f; }
    void SetRingType(RingTypeEnum r)  { ring_ = r; }
  private:
    std::string name_;       ///< Full sugar name
    //std::string resname_;    ///< PDB residue name
    std::string glycamCode_; ///< Glycam residue code
    FormTypeEnum form_;      ///< Sugar form
    ChirTypeEnum chir_;      ///< Sugar chirality
    RingTypeEnum ring_;      ///< Sugar ring type
};

#endif
