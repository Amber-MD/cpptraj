#ifndef INC_MASKARRAY_H
#define INC_MASKARRAY_H
#include <vector>
// Forward declarations
class Topology;
class AtomMask;
namespace Cpptraj {

/// Hold array of atom masks, intended for 'bymol' or 'byres' type selection
class MaskArray {
    typedef std::vector<AtomMask> Marray;
  public:
    MaskArray() {}

    enum SelectType { BY_ATOM = 0, BY_RESIDUE, BY_MOLECULE };

    typedef Marray::const_iterator const_iterator;

    const_iterator begin() const { return masks_.begin(); }
    const_iterator end()   const { return masks_.end();   }

    void SetType(SelectType typeIn) { type_ = typeIn; }
    /// Given an already set up mask, set up internal masks
    int SetupMasks(AtomMask const&, Topology const&);

    unsigned int Nmasks() const { return masks_.size(); }
  private:
    Marray masks_;    ///< Array containing mask selections by residue/molecule
    SelectType type_; ///< The selection type
};

}
#endif
