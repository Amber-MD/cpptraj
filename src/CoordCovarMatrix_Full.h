#ifndef INC_COORDCOVARMATRIX_FULL_H
#define INC_COORDCOVARMATRIX_FULL_H
#include "CoordCovarMatrix.h"
class Frame;
/// Coordinate covanriance full matrix
class CoordCovarMatrix_Full : public CoordCovarMatrix {
  public:
    CoordCovarMatrix_Full();
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    int FinishMatrix();
    // ---------------------------------
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

    Varray vect_1_;
    Varray vect_2_;
    Varray vect2_1_; ///< Array of average diagonal elements squared
    Varray vect2_2_;
    Darray mass1_;
    Darray mass2_;
};
#endif
