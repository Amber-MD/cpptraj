#ifndef INC_COORDCOVARMATRIX_H
#define INC_COORDCOVARMATRIX_H
#include "Matrix.h"
#include <vector>
class Atom;
class AtomMask;
class CpptrajFile;
class DataSet_1D;
class Frame;
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

    /// Print matrix elements to STDOUT for debug with default format
    void DebugPrint(const char*, CpptrajFile&) const;
    /// Print matrix elements to STDOUT for debug with given format
    void DebugPrint(const char*, CpptrajFile&, const char*) const;
  protected:
    typedef Matrix<double> MatType;
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet_1D*> DSarray;

    /// clear internal variables
    virtual void clearMat() = 0;

    /// set mass array
    void set_mass_array(Darray&, std::vector<Atom> const&, AtomMask const&, bool);

    /// \return True if incoming Darray size is divisible by nelt_
    bool has_valid_size(Darray const&) const;

    /// Place coordinates from selected atoms in Frame into an array
    void get_frame_coords(Darray&, Frame const&, AtomMask const&) const;
  //private: // TODO all private

    MatType covarMatrix_;  ///< Coordinate covariance matrix
    unsigned int nframes_; ///< Number of frames added to the matrix
    unsigned int nelt_;    ///< Number of elements per i,j entry
    bool useMass_;         ///< If true use mass weighting
};
#endif
