#ifndef INC_IMAGEROUTINES_H
#define INC_IMAGEROUTINES_H
#include <string>
class Topology;
class Frame;
class AtomMask;
class Matrix_3x3;
class Vec3;
class Box;
#include "ImageTypes.h"
namespace Image {
  class List;
  /// Create an atom pair list by molecule, residue, or atom.
  List* CreateImageList(Topology const&, Mode, std::string const&, bool, bool);
  /// \return Coordinates of center if wrapping molecules into truncated oct. shape.
  Vec3 SetupTruncoct( Frame const&, AtomMask*, bool, bool);
  /// Perform non-orthogonal imaging on given pair list, optionally into trunc. oct. shape.
  void Nonortho(Frame&, bool, Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&,
                bool, List const&);
  /// Perform non-orthogonal imaging on given coordinates, optionally into trunc. oct. shape.
  Vec3 Nonortho(Vec3 const&, bool, bool,
                Matrix_3x3 const&, Matrix_3x3 const&, Vec3 const&, double);
  /// Setup box boundaries based on centering type.
  int SetupOrtho(Box const&, Vec3&, Vec3&, bool);
  /// Perform orthogonal imaging on given pair list using given box boundaries
  void Ortho(Frame&, Vec3 const&, Vec3 const&, Vec3 const&, List const&);
  /// Perform orthogonal imaging on given coordinates using given box boundaries
  Vec3 Ortho(Vec3 const&, Vec3 const&, Vec3 const&, Box const&);
  /// Perform unwrap of non-orthogonal cell using given reference.
  void UnwrapNonortho( Frame&, Frame&, PairType const&, 
                       Matrix_3x3 const&, Matrix_3x3 const&, bool, bool );
  /// Perform unwrap of orthogonal cell using given reference.
  void UnwrapOrtho( Frame&, Frame&, PairType const&, bool, bool );
  /// Wrap selected atom coords from given frame into primary cell, store in result.
  void WrapToCell0(std::vector<double>&, Frame const&, AtomMask const&,
                   Matrix_3x3 const&, Matrix_3x3 const&);
}
#endif
