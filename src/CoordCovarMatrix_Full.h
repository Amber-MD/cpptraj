#ifndef INC_COORDCOVARMATRIX_FULL_H
#define INC_COORDCOVARMATRIX_FULL_H
#include "CoordCovarMatrix.h"
/// Coordinate covanriance full matrix
class CoordCovarMatrix_Full : public CoordCovarMatrix {
  public:
    CoordCovarMatrix_Full();
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    int FinishMatrix();
    // ---------------------------------
    /// Set up matrix for 2 masks
    int SetupMatrix(std::vector<Atom> const&, AtomMask const&,
                    std::vector<Atom> const&, AtomMask const&, bool);
    /// Add selected atoms in Frame to matrix
    void AddFrameToMatrix(Frame const&, AtomMask const&,
                          Frame const&, AtomMask const&);
  private:
    /// Clear the matrix
    void clearMat();
    /// Set up the covariance matrix for selected atoms
    int setupMat(std::vector<Atom> const&, AtomMask const&);
    /// Add elements to the matrix
    void AddToMatrix(Darray const&, Darray const&);

    Darray vect_1_; ///< Hold average of elements in mask1
    Darray vect_2_; ///< Hold average of elements in mask2
    Darray mass1_; ///< Hold masses of atoms in mask1
    Darray mass2_; ///< Hold masses of atoms in mask2
};
#endif
