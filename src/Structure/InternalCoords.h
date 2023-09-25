#ifndef INC_STRUCTURE_INTERNALCOORDS_H
#define INC_STRUCTURE_INTERNALCOORDS_H
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
    /// CONSTRUCTOR - Take pointer to XYZ coords (for first seed atom)
    InternalCoords(const double*);
    /// COPY CONSTRUCTOR
    InternalCoords(InternalCoords const&);
    /// ASSIGNMENT
    InternalCoords& operator=(InternalCoords const&);

    enum IdxType { DISTANCE = 0, ANGLE, TORSION };

    static const int NO_ATOM;
  private:
    int idx_[3];    ///< Atom index for distance, angle, torsion
    double val_[3]; ///< Value for distance, angle, torsion
    double xyz_[3]; ///< XYZ coordinates
    bool isSet_;    ///< True if xyz coords are set
};
}
}
#endif
