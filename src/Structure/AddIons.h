#ifndef INC_CPPTRAJ_STRUCTURE_ADDIONS_H
#define INC_CPPTRAJ_STRUCTURE_ADDIONS_H
#include "../Random.h"
#include <string>
class DataSet_Coords;
class DataSetList;
class Frame;
class Topology;
class Vec3;
namespace Cpptraj {
namespace Parm {
class ParameterSet;
}
namespace Structure {
/// Used to add ions to a system
class AddIons {
  public:
    AddIons();
    /// Init - unit, ion1, Nion1, ion2, Nion2, separation, random seed, debug
    int InitAddIons(std::string const&, int, std::string const&, int, double, int, int);
    /// Print info to stdout
    void PrintAddIonsInfo() const;
    /// Add ions randomly
    int AddIonsRand(Topology&, Frame&, DataSetList const&, Cpptraj::Parm::ParameterSet const&) const;

    /// \return True if at least 1 ion has been specified.
    bool IsSetup() const { return !ion1name_.empty(); }
  private:
    typedef std::vector<Vec3> Varray;

    /// \return Ion unit specified by given name
    DataSet_Coords* GetIonUnit(std::string const&, DataSetList const&) const;
    /// Place ion within topology
    int place_ion(int&, int&, Varray&, DataSet_Coords*, Topology&, Frame&, std::vector<int>&, double, int) const;

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
