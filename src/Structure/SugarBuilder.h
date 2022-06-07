#ifndef INC_STRUCTURE_SUGARBUILDER_H
#define INC_STRUCTURE_SUGARBUILDER_H
#include "SugarToken.h" // FormTypeEnum
#include "../AtomMap.h"
#include <map>
class CpptrajFile;
class NameType;
namespace Cpptraj {
namespace Structure {
// Forward declares
class Sugar;
class ResStatArray;
/// Class for preparing sugars in a topology
class SugarBuilder {
  public:
    /// Type for array of Sugar
    typedef std::vector<Sugar> Array;
    /// CONSTRUCTOR - Take debug level
    SugarBuilder(int);
  private:
    typedef std::vector<int> Iarray;

    /// Set a reduced PDB res to glycam map when dat file not found.
    void SetGlycamPdbResMap();
    /// Load PDB res to glycam map from dat file
    int LoadGlycamPdbResMap(std::string const&);

    /// Find remaining non-ring carbons in chain starting from ring end atom.
    int FindRemainingChainCarbons(Iarray&, int, Topology const&, int, Iarray const&) const;

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
    int ChangePdbAtomNamesToGlycam(std::string const&, Residue const&,
                                   Topology&, Cpptraj::Structure::SugarToken::FormTypeEnum) const;

    /// Determine form/chirality for furanose
    int DetermineUpOrDown(SugarToken&, Sugar const&, Topology const&, Frame const&) const;
    /// Determine form/chirliaty for pyranose 
    int DetermineAnomericForm(SugarToken&, Sugar&, Topology const&, Frame const&) const;

    /// Determine linkages for the sugar
    std::string DetermineSugarLinkages(Sugar const&, CharMask const&, Topology&,
                                       Cpptraj::Structure::ResStatArray&,
                                       std::set<BondType>&,
                                       std::set<BondType>&) const;
    /// Create a residue mask string for selecting Glycam-named sugar residues.
    std::string GenGlycamResMaskString() const;

    /// Try to identify sugar name, form, and linkages
    int IdentifySugar(Sugar&, Topology&, Frame const&, CharMask const&,
                      Cpptraj::Structure::ResStatArray&,
                      std::set<BondType>&,
                      std::set<BondType>&);
    /// Try to find missing linkages to anomeric carbon in sugar.
    int FindSugarC1Linkages(int, int, Topology&, Frame const&, NameType const&) const;
    
    /// ID sugar rings, find missing C1 links, split off functional groups
    int FixSugarsStructure(std::string const&, Topology&, Frame&,
                           bool, bool, NameType const&);
    /// Identify sugars, do renaming, remove bonds, generate leap input
    int PrepareSugars(std::string const&, std::string const&, bool,
                      ResStatArray&,
                      Topology&, Frame const&, CpptrajFile*);

    typedef std::pair<NameType, SugarToken> PairType;
    typedef std::map<NameType, SugarToken> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to sugar information tokens

    typedef std::pair<NameType, NameType> NamePairType;
    typedef std::map<NameType, NameType> NameMapType;

    /// Hold maps of pdb atom names to glycam atom names; multiple residues may share a map
    std::vector<NameMapType> pdb_glycam_name_maps_;
    /// Hold maps of pdb atom names to glycam atom names (res in alpha form)
    std::vector<NameMapType> pdb_glycam_name_maps_A_;
    /// Hold maps of pdb atom names to glycam atom names (res in beta form)
    std::vector<NameMapType> pdb_glycam_name_maps_B_;

    typedef std::pair<std::string, int> ResIdxPairType;
    typedef std::map<std::string, int> ResIdxMapType;
    /// Map glycam residue chars to pdb-glycam atom name map index (pdb_glycam_name_maps_)
    ResIdxMapType glycam_res_idx_map_;

    /// Map pdb residue names to glycam linkage residue names
    NameMapType pdb_glycam_linkageRes_map_;

    Array Sugars_;          ///< Array of found sugars

    bool hasGlycam_;        ///< If true, assume sugars already have glycam names
    bool useSugarName_;     ///< If true, base form/chirality on name instead of geometry
    AtomMap myMap_;         ///< Used to determine unique atoms for chirality
    int debug_;
};
}
}
#endif
