#ifndef INC_VIEWRST_H
#define INC_VIEWRST_H
#include "Topology.h"
/// Used for creating mol2 that can be used to visualize restraints.
class ViewRst {
  public:
    enum NoeType {NOE_STRONG=0, NOE_MEDIUM, NOE_WEAK, NOE_VERYWEAK};
    enum OutputType { ALL=0, BY_STRENGTH };
    ViewRst() : outType_(ALL) {}
    /// Initialize with given topology, frame, and output type.
    int Init(Topology const&, Frame const&, OutputType);
    /// Add restraint, default strong (for ALL)
    void AddRst(int, int, NoeType);
    /// Add restraint of given strength (for BY_STRENGTH)
    void AddRst(int i, int j) { AddRst(i, j, NOE_STRONG); }
    /// Write output mol2 file(s) using given file name.
    int WriteRstMol2(std::string const&);
  private:
    std::vector< Topology > Pseudo_;      ///< Pseudo topology(s) containing rst bonds
    Frame coords_;                        ///< Output coordinates
    OutputType outType_;                  ///< Output type.
};
#endif
