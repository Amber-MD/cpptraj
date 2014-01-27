#ifndef INC_SYMMETRICRMSDCALC_H
#define INC_SYMMETRICRMSDCALC_H
#include "Topology.h"
#include "Hungarian.h"
/// Class for performing symmetric RMSD calculations.
class SymmetricRmsdCalc {
  public:
    SymmetricRmsdCalc();
    int FindSymmetricAtoms(Topology const&, AtomMask const&);
    double SymmRMSD(Frame const&, AtomMask const&, Frame const&, Frame const&,
                    Matrix_3x3&, Vec3&, Vec3 const&, bool, bool);
    const Frame* RemapFrame() const { return &remapFrame_; }
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> AtomIndexArray;
    /// Array of groups of potentially symmetric atoms
    AtomIndexArray SymmetricAtomIndices_;
    int debug_;
    Hungarian cost_matrix_;
    Iarray AMap_;           ///< AMap_[ref] = tgt
    Frame remapFrame_;      ///< Target frame re-mapped for symmetry
    Frame tgtFrame_;        ///< Selected atoms from target frame.
};
#endif
