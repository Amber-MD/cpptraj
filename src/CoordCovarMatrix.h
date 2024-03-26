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
    /// Clear the matrix
    void Clear();

    /// Print matrix elements to STDOUT for debug
    void DebugPrint(const char*, CpptrajFile&) const;
  protected:
    typedef Matrix<double> MatType;
    typedef std::vector<double> Darray;
    typedef std::vector<Vec3> Varray;

    /// clear internal variables
    virtual void clearMat() = 0;

    /// set mass array
    void set_mass_array(Darray&, std::vector<Atom> const&, AtomMask const&, bool);

  //private: // TODO all private

    MatType covarMatrix_;  ///< Coordinate covariance matrix
    unsigned int nframes_; ///< Number of frames added to the matrix
    bool useMass_;         ///< If true use mass weighting
};
#endif
