#ifndef INC_STRUCTURE_INTERNALCOORDS_H
#define INC_STRUCTURE_INTERNALCOORDS_H
class Vec3;
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

    //enum ValType { DISTANCE = 0, ANGLE, TORSION };
    //enum AtType  { AT_J = 0,     AT_K,  AT_L };

    static const int NO_ATOM;

    /// Zero XYZ coords (seed atom 0)
    void ZeroXYZ();
    /// Set XYZ coords
    void SetXYZ( Vec3 const&);
    /// \return Specifed value
    //double Val(ValType v) const { return val_[(int)v]; }
    /// \return Specified atom index
    //int Idx(AtType a) const { return idx_[(int)a]; }
    double Dist() const { return val_[0]; }
    double Theta() const { return val_[1]; }
    double Phi() const { return val_[2]; }

    int AtJ() const { return idx_[0]; }
    int AtK() const { return idx_[1]; }
    int AtL() const { return idx_[2]; }

    /// Set index of atom
    //void SetIdx(AtType a) { idx_[(int)a] = idx; }
  private:
    int idx_[3];    ///< Atom index for distance, angle, torsion
    double val_[3]; ///< Value for distance, angle, torsion
    double xyz_[3]; ///< XYZ coordinates
    bool isSet_;    ///< True if xyz coords are set
};
}
}
#endif
