#ifndef INC_COORDCOVARMATRIX_H
#define INC_COORDCOVARMATRIX_H
#include "Matrix.h"
#include "Vec3.h"
#include <vector>
class Atom;
class AtomMask;
class CpptrajFile;
class Frame;
/// Coordinate covariance matrix
class CoordCovarMatrix {
  public:
    /// CONSTRUCTOR
    CoordCovarMatrix();
    /// Clear the matrix
    void Clear();
    /// Set up the covariance matrix for selected atoms
    int SetupMatrix(std::vector<Atom> const&, AtomMask const&, bool);
    /// Add selected atoms in Frame to matrix
    void AddFrameToMatrix(Frame const&, AtomMask const&);
    /// Add Frame to matrix
    void AddFrameToMatrix(Frame const&);
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    int FinishMatrix();

    /// Print matrix elements to STDOUT for debug
    void DebugPrint(const char*, CpptrajFile&) const;
  private:
    typedef Matrix<double> MatType;
    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;

    MatType covarMatrix_;  ///< Coordinate covariance matrix
    Varray vect_;          ///< Store average coordinates along the diagonal
    Darray mass_;          ///< Store selected atoms masses
    unsigned int nframes_; ///< Number of frames added to the matrix
    bool useMass_;         ///< If true use mass weighting
};
#endif
