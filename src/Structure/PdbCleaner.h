#ifndef INC_STRUCTURE_PDBCLEANER_H
#define INC_STRUCTURE_PDBCLEANER_H
#include <vector>
#include <string>
class ArgList;
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Used to clean up structures read from a PDB
class PdbCleaner {
    typedef std::vector<int> Iarray;
  public:
    PdbCleaner();
    void SetDebug(int d) { debug_ = d; }
    /// Help text
    static void Help();
    /// Initialize
    int InitPdbCleaner(ArgList&, std::string const&, Iarray const&);
    /// Setup
    int SetupPdbCleaner(Topology const&);
    /// Info
    void PdbCleanerInfo() const;
    /// Modify top/frame
    int ModifyCoords(Topology&, Frame&) const;
    /// \return True if option to remove water is set
    bool RemoveWater() const { return remove_water_; }
  private:
    int debug_;
    bool remove_water_;      ///< If true, remove any water.
    bool remove_h_;          ///< If true, remove hydrogen atoms.
    std::string waterMask_;  ///< Mask expression for selecting water.
    std::string altLocArg_;  ///< Only keep atoms with this alternate location identifier.
    std::string stripMask_;  ///< General mask for removing atoms.
    Iarray resnumsToRemove_; ///< Other residue numbers to remove.
};
}
}
#endif
