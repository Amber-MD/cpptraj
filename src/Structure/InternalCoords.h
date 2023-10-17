#ifndef INC_STRUCTURE_INTERNALCOORDS_H
#define INC_STRUCTURE_INTERNALCOORDS_H
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
  private:
    int ati_;       ///< Atom I index
    int idx_[3];    ///< Atom indices for distance, angle, torsion
    double val_[3]; ///< Values for distance, angle, torsion
};
}
}
#endif
