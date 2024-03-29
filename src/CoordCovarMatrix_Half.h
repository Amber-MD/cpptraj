#ifndef INC_COORDCOVARMATRIX_HALF_H
#define INC_COORDCOVARMATRIX_HALF_H
#include "CoordCovarMatrix.h"
class Frame;
/// Coordinate covariance half (self) matrix
class CoordCovarMatrix_Half : public CoordCovarMatrix {
  public:
    /// CONSTRUCTOR 
    CoordCovarMatrix_Half();
    // ---------------------------------
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    int FinishMatrix();
    // ---------------------------------
    /// Set up half matrix
    int SetupMatrix(std::vector<Atom> const&, AtomMask const&, bool);
    /// Add selected atoms in Frame to matrix
    void AddFrameToMatrix(Frame const&, AtomMask const&);
    /// Add Frame to matrix
    //void AddFrameToMatrix(Frame const&);
  private:
    /// Add elements to the matrix
    void AddToMatrix(Darray const&);
    /// Clear the matrix
    void clearMat();

    Darray vect_;
    Darray mass_;

};
#endif
