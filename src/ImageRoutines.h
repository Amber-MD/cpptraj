#ifndef INC_IMAGEROUTINES_H
#define INC_IMAGEROUTINES_H
#include "Frame.h"
// TODO: Make Frame const
Vec3 SetupImageTruncoct( Frame&, AtomMask*, bool, bool);
void ImageNonortho(Frame&, bool, Vec3 const&, Matrix_3x3, Matrix_3x3,
                   bool, bool, bool, std::vector<int> const&);
Vec3 ImageNonortho(Vec3 const&, bool, bool,
                   const Matrix_3x3&, const Matrix_3x3&, Vec3 const&, double);
void SetupImageOrtho(Frame&, Vec3&, Vec3&, bool);
void ImageOrtho(Frame&, Vec3 const&, Vec3 const&, bool, bool, std::vector<int> const&);
Vec3 ImageOrtho(Vec3 const&, Vec3 const&, Vec3 const&, Vec3 const&);
void UnwrapNonortho( Frame&, Frame&, AtomMask const&, Matrix_3x3 const&, Matrix_3x3 const& );
void UnwrapOrtho( Frame&, Frame&, AtomMask const& );
#endif
