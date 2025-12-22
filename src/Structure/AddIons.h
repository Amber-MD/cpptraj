#ifndef INC_CPPTRAJ_STRUCTURE_ADDIONS_H
#define INC_CPPTRAJ_STRUCTURE_ADDIONS_H
#include "../Random.h"
#include <string>
namespace Cpptraj {
namespace Structure {
/// Used to add ions to a system
class AddIons {
  public:
    AddIons();
    /// Init - unit, ion1, Nion1, ion2, Nion2, separation, random seed, debug
    int InitAddIons(std::string const&, std::string const&, int, std::string const&, int, double, int, int);
  private:
    std::string unitname_; ///< System to add ions to
    std::string ion1name_; ///< Name of first ion
    std::string ion2name_; ///< Name of second ion
    int Nion1_; ///< Number of first ion to add
    int Nion2_; ///< Number of second ion to add
    int debug_;
    double separation_; ///< Min. distance in Angstroms any ion can be from another ion
    Random_Number RNG_; ///< Random number generator
};
}
}
#endif
