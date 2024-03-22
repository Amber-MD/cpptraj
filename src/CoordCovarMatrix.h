#ifndef INC_COORDCOVARMATRIX_H
#define INC_COORDCOVARMATRIX_H
#include "Matrix.h"
#include "Vec3.h"
#include <vector>
class Atom;
class AtomMask;
class CpptrajFile;
/// Coordinate covariance matrix abstract base class
class CoordCovarMatrix {
  public:
    /// CONSTRUCTOR
    CoordCovarMatrix();
    /// DESTRUCTOR - virtual since inherited
    virtual ~CoordCovarMatrix() {}
    // --------------------------------- 
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    virtual int FinishMatrix() = 0;
    // ---------------------------------
    /// Set up the covariance matrix for selected atoms
    int SetupMatrix(std::vector<Atom> const&, AtomMask const&, bool);
    /// Clear the matrix
    void Clear();

    /// Print matrix elements to STDOUT for debug
    void DebugPrint(const char*, CpptrajFile&) const;
  protected:
    /// clear internal variables
    virtual void clearMat() = 0;
    /// set internal variables
    virtual int setupMat(std::vector<Atom> const&, AtomMask const&) = 0;
  //private: // TODO all private
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
