#ifndef INC_COORDCOVARMATRIX_H
#define INC_COORDCOVARMATRIX_H
#include "Matrix.h"
#include "Vec3.h"
#include <vector>
class AtomMask;
class Frame;
/// Coordinate covariance matrix
class CoordCovarMatrix {
  public:
    /// CONSTRUCTOR
    CoordCovarMatrix();
    /// Clear the matrix
    void Clear();
    /// Add Frame to matrix
    void AddFrame(Frame const&, AtomMask const&);
  private:
    typedef Matrix<double> MatType;
    typedef std::vector<Vec3> Varray;

    MatType covarMatrix_;  ///< Coordinate covariance matrix
    Varray vect_;          ///< Store average coordinates along the diagonal
    unsigned int nframes_; ///< Number of frames added to the matrix
};
#endif
