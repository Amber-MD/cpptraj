#ifndef INC_SYMMETRICRMSDCALC_H
#define INC_SYMMETRICRMSDCALC_H
#include "Topology.h"
#include "Hungarian.h"
/// Class for performing symmetric RMSD calculations.
class SymmetricRmsdCalc {
  public:
    SymmetricRmsdCalc();
    /// Set target mask string, fit, and mass options.
    int InitSymmRMSD(std::string const&, bool, bool, int);
    /// Setup target mask, find symmetric atoms.
    int FindSymmetricAtoms(Topology const&);
    double SymmRMSD(Frame const&, Frame const&, Frame const&,
                    Vec3 const&, Matrix_3x3&, Vec3&);
    AtomMask const& TgtMask() const { return tgtMask_;     }
    const Frame* RemapFrame() const { return &remapFrame_; }
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> AtomIndexArray;
    /// Array of groups of potentially symmetric atoms
    AtomIndexArray SymmetricAtomIndices_;
    int debug_;
    Hungarian cost_matrix_; ///< Hungarian algorithm cost matrix.
    Iarray AMap_;           ///< AMap_[ref] = tgt
    Frame remapFrame_;      ///< Target frame re-mapped for symmetry
    Frame selectedTgt_;     ///< Selected atoms from target frame.
    AtomMask tgtMask_;      ///< Mask selecting atoms in target for RMSD calc.
    bool fit_;              ///< If true, perform RMS best-fit.
    bool useMass_;          ///< If true, mass-weight calc.
};
#endif
