#ifndef INC_PUCKER_PUCKERMASK_H
#define INC_PUCKER_PUCKERMASK_H
#include <vector>
#include <string>
class Topology;
namespace Cpptraj {
namespace Pucker {
/// Used to define found puckers in Topology
class PuckerMask {
  public:
    PuckerMask();
    /// CONSTRUCTOR - residue number, name (DataSet aspect), res # (DataSet index)
    PuckerMask(int, std::string const&, std::vector<int> const&);
    /// \return true if mask is empty
    bool None() const { return atoms_.empty(); }
    /// \return string based on atoms in the mask
    std::string PuckerMaskString(Topology const&) const;
    /// \return Residue number pucker belongs to
    int ResNum() const { return resnum_; }
    /// \return Name of pucker
    std::string const& Name() const { return aspect_; }
  private:
    std::vector<int> atoms_; ///< Hold atom indices defining pucker
    std::string aspect_;     ///< DataSet aspect
    int resnum_;             ///< Residue number (DataSet index)
};

}
}
#endif
