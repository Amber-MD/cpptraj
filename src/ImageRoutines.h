#ifndef INC_IMAGEROUTINES_H
#define INC_IMAGEROUTINES_H
#include <string>
#include <vector>
#include "ImageTypes.h"
class Topology;
class Frame;
class AtomMask;
class Matrix_3x3;
class Vec3;
class Box;
class Unit;

namespace Image {
  class List;
  /// \return empty image list for given mode, optionally using mass and centering
  List* CreateImageList(Mode, bool, bool);
  /// \return an image list for given Top, mode, mask, optionally using mass and centering
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
  void UnwrapNonortho( Frame&, Frame&, List const&, Unit const&,
                       Matrix_3x3 const&, Matrix_3x3 const& );
  /// Perform unwrap of orthogonal cell using given reference.
  void UnwrapOrtho( Frame&, Frame&, List const&, Unit const& );

  /// Wrap selected atom coords from given frame into primary cell, store in result.
  void WrapToCell0(std::vector<double>&, Frame const&, AtomMask const&,
                   Matrix_3x3 const&, Matrix_3x3 const&);
}
#endif
