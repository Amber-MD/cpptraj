#ifndef INC_SYMMETRICRMSDCALC_H
#define INC_SYMMETRICRMSDCALC_H
#include "Topology.h"
#include "Hungarian.h"
#include "AtomMap.h"
/// Class for performing symmetry-corrected RMSD calculations.
class SymmetricRmsdCalc {
  public:
    SymmetricRmsdCalc();
    /// Set target mask string, fit, and mass options.
    int InitSymmRMSD(std::string const&, bool, bool, int);
    /// Setup target mask, find symmetric atoms.
    int SetupSymmRMSD(Topology const&);
    /// Calculate symmetry-corrected RMSD using both original and pre-centered references.
    double SymmRMSD(Frame const&, Frame const&, Frame const&, Vec3 const&);
    AtomMask const& TgtMask()     const { return tgtMask_;     }
    const Frame* RemapFrame()     const { return &remapFrame_; }
    bool Fit()                    const { return fit_;         }
    bool UseMass()                const { return useMass_;     }
    Matrix_3x3 const& RotMatrix() const { return rotMatrix_;   }
    Vec3 const& TgtTrans()        const { return tgtTrans_;    }
  private:
    enum atomStatusType { UNSELECTED = 0, NONSYMM, SYMM };
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> AtomIndexArray;
    
    void FindSymmetricAtoms(int, AtomMap const&, std::string const&, Iarray&, Iarray&) const;

    /// Array of groups of potentially symmetric atoms
    AtomIndexArray SymmetricAtomIndices_;
    int debug_;
    Hungarian cost_matrix_; ///< Hungarian algorithm cost matrix.
    Iarray AMap_;           ///< AMap_[ref] = tgt
    Frame remapFrame_;      ///< Target frame re-mapped for symmetry
    Frame selectedTgt_;     ///< Selected atoms from target frame.
    AtomMask tgtMask_;      ///< Mask selecting atoms in target for RMSD calc.
    Matrix_3x3 rotMatrix_;  ///< Hold best-fit rotation matrix for target.
    Vec3 tgtTrans_;         ///< Hold translation of target to origin.
    bool fit_;              ///< If true, perform RMS best-fit.
    bool useMass_;          ///< If true, mass-weight calc.
};
#endif
