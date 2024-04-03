#ifndef INC_COORDCOVARMATRIX_HALF_H
#define INC_COORDCOVARMATRIX_HALF_H
#include "CoordCovarMatrix.h"
class DataSet_1D;
/// Coordinate covariance half (self) matrix
class CoordCovarMatrix_Half : public CoordCovarMatrix {
  public:
    /// CONSTRUCTOR 
    CoordCovarMatrix_Half();
    // ---------------------------------
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    int FinishMatrix();
    // ---------------------------------
    /// Set up half matrix for coordinates
    int SetupMatrix(std::vector<Atom> const&, AtomMask const&, bool);
    /// Set up half matrix for data sets
    int SetupMatrix(std::vector<DataSet_1D*> const&);
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
