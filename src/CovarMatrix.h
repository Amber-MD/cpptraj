#ifndef INC_COVARMATRIX_H
#define INC_COVARMATRIX_H
#include "Matrix.h"
#include "Vec3.h"
#include <vector>
class Atom;
class AtomMask;
class CpptrajFile;
class DataSet_2D;
/// Coordinate covariance matrix abstract base class
class CovarMatrix {
  public:
    /// CONSTRUCTOR - Number of elements
    CovarMatrix(unsigned int);
    /// DESTRUCTOR - virtual since inherited
    virtual ~CovarMatrix() {}
    // --------------------------------- 
    /// Finish calculating the matrix (normalize, calc <rirj> - <ri><rj>)
    virtual int FinishMatrix() = 0;
    // ---------------------------------
    /// Clear the matrix
    void Clear();

    /// Print matrix elements to STDOUT for debug
    void DebugPrint(const char*, CpptrajFile&) const;
  protected:
    typedef DataSet_2D* MatType;
    typedef std::vector<double> Darray;

    /// clear internal variables
    virtual void clearMat() = 0;

    /// set mass array
    void set_mass_array(Darray&, std::vector<Atom> const&, AtomMask const&, bool);

  //private: // TODO all private

    MatType covarMatrix_;   ///< Covariance matrix
    unsigned int nelt_;     ///< Number of elements per atom
    unsigned int nframes_;  ///< Number of frames added to the matrix
    bool useMass_;          ///< If true use mass weighting
};
#endif
