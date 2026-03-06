#ifndef INC_STRUCTURE_METALCENTERFINDER_H
#define INC_STRUCTURE_METALCENTERFINDER_H
#include "../AtomMask.h"
#include <vector>
#include <map>
class ArgList;
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Used for finding and preparing metal centers
class MetalCenterFinder {
  public:
    /// CONSTRUCTOR
    MetalCenterFinder();
    /// Init with arguments
    int InitMetalCenters(ArgList&, int);
    /// Find metal centers
    int FindMetalCenters(Topology const&, Frame const&);
    /// Print Info to stdout
    void PrintMetalCenterInfo() const;
    /// Print found metal centers to stdout
    void PrintMetalCenters(Topology const&) const;
  private:
    typedef std::vector<int> Iarray;
    typedef std::pair<int, Iarray> MCpair;
    /// Map metal center index to coordinating atom indices
    typedef std::map<int, Iarray> MCmap;

    MCmap metalCenters_;     ///< Hold metal centers and coordinating atoms
    AtomMask metalMask_;     ///< Mask containing potential metal centers
    AtomMask coordAtomMask_; ///< Mask containing potential coordinating atoms
    double dcut2_;           ///< Distance cutoff in Ang^2
    int debug_;
};
}
}
#endif
