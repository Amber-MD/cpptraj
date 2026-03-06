#ifndef INC_PARM_ASSIGNPARAMS_H
#define INC_PARM_ASSIGNPARAMS_H
#include <vector>
#include <string>
#ifdef TIMER
# include "../Timer.h"
#endif
class AngleArray;
class AngleParmArray;
class AngleParmType;
class AngleType;
class Atom;
class AtomType;
class BondArray;
class BondParmArray;
class BondParmType;
class BondType;
class CmapArray;
class CmapGridArray;
class CmapParmHolder;
class DihedralArray;
class DihedralParmArray;
class DihedralType;
class HB_ParmType;
class NonbondType;
class Topology;
namespace Cpptraj {
namespace Parm {
// Parm forward declares
template<typename Type> class ParmHolder;
class DihedralParmHolder;
class ImproperParmHolder;
class ParameterSet;
/// Used to assign parameters to a Topology
class AssignParams {
  public:
    AssignParams();
    /// Set Debug level
    void SetDebug(int);
    /// Set verbosity
    void SetVerbose(int);
    /// Allow flexible water
    void SetFlexibleWater(bool);
    /// Replace existing parameters with those from given set
#   ifdef TIMER
    int AssignParameters(Topology&, ParameterSet const&, int&);
#   else
    int AssignParameters(Topology&, ParameterSet const&, int&) const;
#   endif
    /// Update existing parameters with given parameter set
    int UpdateParameters(Topology&, ParameterSet const&)
#   ifdef TIMER
      ;
#   else
      const;
#   endif
    // Assign nonbond parameters to topology. Used during Topology::AppendTop()
    void AssignNonbondParams(Topology&,
                             ParmHolder<AtomType> const&,
                             ParmHolder<NonbondType> const&, ParmHolder<NonbondType> const&,
                             ParmHolder<double> const&, ParmHolder<HB_ParmType> const&) const;
    void WriteAssignTiming(int, double) const;
  private:
    typedef std::vector<Atom> AtArray;

    int AssignAtomTypeParm(AtArray&, ParmHolder<AtomType> const&) const;
    int AssignBondParm(Topology const&, ParmHolder<BondParmType> const&,
                        BondArray&, BondParmArray&, const char*) const;
    void AssignBondParams(Topology&, ParmHolder<BondParmType> const&) const;
    void AssignUBParams(Topology&, ParmHolder<BondParmType> const&) const;

    AngleArray AssignAngleParm(Topology const&, ParmHolder<AngleParmType> const&,
                               AngleArray const&, AngleParmArray&, int&) const;
    void AssignAngleParams(Topology&, ParmHolder<AngleParmType> const&) const;

    void warn_improper_reorder(AtArray const&, DihedralType const&, DihedralType const&) const;
    void AssignImproperParm(Topology&, ImproperParmHolder const&,
                            DihedralArray&, DihedralParmArray&) const ;
    void AssignImproperParams(Topology&, ImproperParmHolder const&) const;

    DihedralArray AssignDihedralParm(Topology&,
                                     DihedralParmHolder const&, ImproperParmHolder const&,
                                     ParmHolder<AtomType> const&, DihedralArray const&,
                                     bool, int&)  // TODO make the bool a class var?
#   ifdef TIMER
      ;
#   else
      const;
#   endif
    DihedralArray get_unique_dihedrals(DihedralArray const&) const;
    void AssignDihedralParams(Topology&, DihedralParmHolder const&, ImproperParmHolder const&,
                              ParmHolder<AtomType> const&)
#   ifdef TIMER
      ;
#   else
      const;
#   endif

    static inline int cmap_anames_match(DihedralType const&, AtArray const&, std::vector<std::string> const&);
    int remap_cmap_indices(Topology const&, std::vector<int>&, CmapGridArray&, CmapArray&, CmapParmHolder const&) const;
    int AssignCmapParams(Topology const&, CmapArray&, CmapParmHolder const&, CmapGridArray&) const;
    int AssignCmapParams(Topology const&, DihedralArray const&, CmapParmHolder const&,
                         CmapGridArray&, CmapArray&) const;
    static void AddToBondArrays(Topology&, BondType const&);
    static void AddToAngleArrays(Topology&, AngleType const&);
    static void AddToDihedralArrays(Topology&, DihedralType const&);

    int debug_;
    int verbose_;
    bool deleteExtraPointAngles_; ///< If true, remove angles/torsions containing extra points.
    bool flexibleWater_;          ///< If true, allow H-O-H angle for water
#   ifdef TIMER
    Timer t_total_;
    Timer t_bonds_;
    Timer t_angles_;
    Timer t_dihedrals_;
    Timer t_cmaps_;
    Timer t_multi_;
    Timer t_dih_imp_;
    Timer t_dih_dih_;
    Timer t_dih_dih_getnew_;
    Timer t_UB_;
    Timer t_improper_;
    Timer t_atype_;
    Timer t_nonbond_;
#   endif
};
}
}
#endif
