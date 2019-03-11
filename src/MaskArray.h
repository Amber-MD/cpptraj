#ifndef INC_MASKARRAY_H
#define INC_MASKARRAY_H
#include <vector>
#include <string>
#include "AtomMask.h"
// Forward declarations
class Topology;
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

    int SetMaskExpression(std::string const&, SelectType);
    int SetupMasks(Topology const&);

    unsigned int Nmasks() const { return masks_.size(); }
  private:
    AtomMask mainMask_; ///< The overall atom mask
    Marray masks_;
    SelectType type_;   ///< The selection type
};

}
#endif
