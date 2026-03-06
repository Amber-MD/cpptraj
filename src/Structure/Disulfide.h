#ifndef INC_STRUCTURE_DISULFIDE_H
#define INC_STRUCTURE_DISULFIDE_H
#include <string>
#include <vector>
class ArgList;
class BondType;
class CpptrajFile;
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
class ResStatArray;
/// Used to search for disulfide bonds and rename residues appropriately
class Disulfide {
  public:
    /// Constructor
    Disulfide();
    /// Keywords recognized by InitDisulfide
    static const char* keywords_;
    /// Will indicate whether we want to add newly detected disulfide bonds to Topology
    enum AddType { ADD_BONDS = 0, NO_ADD_BONDS };
    /// Init from args
    int InitDisulfide(ArgList&, AddType, int);
    /// Print info to stdout
    //void DisulfideInfo() const;
    /// Search for disulfide bonds
    int SearchForDisulfides(ResStatArray&, Topology&, Frame const&, std::vector<BondType>&) const;
  private:
    /// Search for disulfide bonds
    static int searchForDisulfides(ResStatArray&,
                        double, std::string const&, std::string const&, bool,
                        Topology&, Frame const&, std::vector<BondType>&, AddType);

    double disulfidecut_; ///< Distance cutoff for disulfide bonds in Ang.
    std::string newcysnamestr_; ///< Name of cysteine residues involved in disulfide bonds
    std::string cysmaskstr_;    ///< Cysteine mask string
    bool searchForNewDisulfides_; ///< If true, look for new disulfides in addition to existing bonds
    AddType addNewBonds_; ///< Indicate whether newly detected disulfides should be added to Topology
};
}
}
#endif
