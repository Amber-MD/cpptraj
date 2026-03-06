#ifndef INC_PARM_GETPARAMS_H
#define INC_PARM_GETPARAMS_H
#include <vector>
class AngleArray;
class AngleParmArray;
class AngleParmType;
class Atom;
class AtomType;
class BondArray;
class BondParmArray;
class BondParmType;
class CmapArray;
class CmapGridArray;
class CmapParmHolder;
class DihedralArray;
class DihedralParmArray;
class HB_ParmType;
class NonbondParmType;
class NonbondType;
class Residue;
class Topology;
namespace Cpptraj {
namespace Parm {
template<typename Type> class ParmHolder;
class DihedralParmHolder;
class ImproperParmHolder;
class ParameterSet;
/// Used to get parameters from a Topology
class GetParams {
  public:
    /// CONSTRUCTOR
    GetParams();
    /// Set debug level
    void SetDebug(int);
    /// \return ParameterSet containing parameters from given topology.
    ParameterSet GetParameters(Topology const&) const;
    /// Get nonbonded parameters (used in AppendTop)
    void GetLJAtomTypes(ParmHolder<AtomType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<double>&,
                        ParmHolder<HB_ParmType>&,
                        std::vector<Atom> const&,
                        NonbondParmType const&) const;
    /// \return number of unique atom types
    unsigned int NuniqueAtomTypes(Topology const&) const;
  private:
    static inline void GetBondParams(ParmHolder<BondParmType>&, std::vector<Atom> const&,
                                     BondArray const&, BondParmArray const&);
    static inline void GetAngleParams(ParmHolder<AngleParmType>&, std::vector<Atom> const& atoms,
                                      AngleArray const&, AngleParmArray const&);
    static inline void GetImproperParams(ImproperParmHolder&, std::vector<Atom> const&,
                                         DihedralArray const&, DihedralParmArray const&);
    static inline void GetDihedralParams(DihedralParmHolder&, ImproperParmHolder&,
                                         std::vector<Atom> const&,
                                         DihedralArray const&, DihedralParmArray const&);
    inline int GetCmapParams(CmapParmHolder&, CmapArray const&, CmapGridArray const&,
                                    std::vector<Atom> const&, std::vector<Residue> const&) const;

    int debug_;
};
}
}
#endif
