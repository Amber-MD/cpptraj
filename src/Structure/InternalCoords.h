#ifndef INC_STRUCTURE_INTERNALCOORDS_H
#define INC_STRUCTURE_INTERNALCOORDS_H
#include <vector>
class Topology;
namespace Cpptraj {
namespace Structure {
/// Hold internal coordinates for an atom
/** Given that this is atom i connected to atoms j, k, and l:
  *   k - j
  *  /     \
  * l       i
  * hold the following:
  *   distance j - i
  *   angle    k - j - i
  *   torsion  l - k - j - i
  */
class InternalCoords {
  public:
    /// CONSTRUCTOR
    InternalCoords();
    /// CONSTRUCTOR - Take atoms i, j, k, l and values dist, theta, phi
    InternalCoords(int, int, int, int, double, double, double);
    /// COPY CONSTRUCTOR
    InternalCoords(InternalCoords const&);
    /// ASSIGNMENT
    InternalCoords& operator=(InternalCoords const&);
    /// Less than another IC using AI < AL < AK < AJ
    bool operator<(InternalCoords const& rhs) const {
      if (ati_ == rhs.ati_) {
        if (idx_[2] == rhs.idx_[2]) { // Atom L
          if (idx_[1] == rhs.idx_[1]) { // Atom K
            if (idx_[0] == rhs.idx_[0]) { // Atom J
              return false; // Equal
            } else {
              return (idx_[0] < rhs.idx_[0]);
            }
          } else {
            return (idx_[1] < rhs.idx_[1]);
          }
        } else {
          return (idx_[2] < rhs.idx_[2]);
        }
      } else {
        return (ati_ < rhs.ati_);
      }
    }
    /// \return True if all indices are equal
    bool operator==(InternalCoords const& rhs) const {
      return ( ati_ == rhs.ati_ &&
               idx_[0] == rhs.idx_[0] &&
               idx_[1] == rhs.idx_[1] &&
               idx_[2] == rhs.idx_[2] );
    }
    /// \return True if any index is not equal
    bool operator!=(InternalCoords const& rhs) const {
      return ( ati_ != rhs.ati_ ||
               idx_[0] != rhs.idx_[0] ||
               idx_[1] != rhs.idx_[1] ||
               idx_[2] != rhs.idx_[2] );
    }
    /// Indictaes no atom set
    static const int NO_ATOM;

    double Dist() const { return val_[0]; }
    double Theta() const { return val_[1]; }
    double Phi() const { return val_[2]; }

    int AtI() const { return ati_; }
    int AtJ() const { return idx_[0]; }
    int AtK() const { return idx_[1]; }
    int AtL() const { return idx_[2]; }

    static unsigned int sizeInBytes() { return (4*sizeof(int)) + (3*sizeof(double)); }

    void printIC(Topology const&) const;
    /// Remap internal indices according to given map
    //void RemapIndices(std::vector<int> const&);
    /// Offset internal indices by given offset
    void OffsetIndices(int);
    /// Set the value of Phi
    void SetPhi(double p) { val_[2] = p; }
    /// Set the value of theta
    void SetTheta(double t) { val_[1] = t; }
    /// Set the value of distance
    void SetDist(double d) { val_[0] = d; }
  private:
    int ati_;       ///< Atom I index
    int idx_[3];    ///< Atom indices for distance, angle, torsion
    double val_[3]; ///< Values for distance, angle, torsion
};
}
}
#endif
