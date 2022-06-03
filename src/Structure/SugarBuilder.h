#ifndef INC_STRUCTURE_SUGARBUILDER_H
#define INC_STRUCTURE_SUGARBUILDER_H
#include "SugarToken.h" // FormTypeEnum
#include "../AtomMap.h"
#include <map>
namespace Cpptraj {
namespace Structure {
// Forward declares
class Sugar;
class ResStatArray;
/// Class for preparing sugars in a topology
class SugarBuilder {
  public:
    SugarBuilder();
  private:
    typedef std::vector<int> Iarray;

    /// Set a reduced PDB res to glycam map when dat file not found.
    void SetGlycamPdbResMap();
    /// Load PDB res to glycam map from dat file
    int LoadGlycamPdbResMap(std::string const&);

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
                                       CpptrajFile*, std::set<BondType>&) const;
    /// Create a residue mask string for selecting Glycam-named sugar residues.
    std::string GenGlycamResMaskString() const;
    /// Try to identify sugar name, form, and linkages
    int IdentifySugar(Sugar&, Topology&, Frame const&, CharMask const&, CpptrajFile*, std::set<BondType>&);
    /// Try to find missing linkages to anomeric carbon in sugar.
    int FindSugarC1Linkages(int, int, Topology&, Frame const&) const;
    
    /// Attempt to fix any issues with sugars
    int FixSugarsStructure(std::vector<Sugar>&, std::string const&, Topology&, Frame&,
                           bool, bool) const;

    int PrepareSugars(std::string const&, std::vector<Sugar>&, Topology&, Frame const&, CpptrajFile*);

    typedef std::pair<NameType, SugarToken> PairType;
    typedef std::map<NameType, SugarToken> MapType;
    MapType pdb_to_glycam_; ///< Map PDB residue names to sugar information tokens
 

    bool hasGlycam_;        ///< If true, assume sugars already have glycam names
    bool useSugarName_;     ///< If true, base form/chirality on name instead of geometry
    AtomMap myMap_;         ///< Used to determine unique atoms for chirality

#endif
