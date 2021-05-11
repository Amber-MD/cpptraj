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

    PuckerMask(std::vector<int> const&);

    bool None() const { return atoms_.empty(); }
    std::string PuckerMaskString(Topology const&) const;
  private:
    std::vector<int> atoms_; ///< Hold atom indices defining pucker
};

}
}
#endif
