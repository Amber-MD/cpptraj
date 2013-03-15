#ifndef INC_IMAGEROUTINES_H
#define INC_IMAGEROUTINES_H
#include "Frame.h"
Vec3 SetupImageTruncoct( Frame const&, AtomMask*, bool, bool);
void ImageNonortho(Frame&, bool, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&,
                   bool, bool, bool, std::vector<int> const&);
Vec3 ImageNonortho(Vec3 const&, bool, bool,
                   Matrix_3x3 const&, Matrix_3x3 const&, Vec3 const&, double);
int SetupImageOrtho(Box const&, Vec3&, Vec3&, bool);
void ImageOrtho(Frame&, Vec3 const&, Vec3 const&, bool, bool, std::vector<int> const&);
Vec3 ImageOrtho(Vec3 const&, Vec3 const&, Vec3 const&, Box const&);
void UnwrapNonortho( Frame&, Frame&, AtomMask const&, Matrix_3x3 const&, Matrix_3x3 const& );
void UnwrapOrtho( Frame&, Frame&, AtomMask const& );
#endif
