#ifndef INC_BOXTYPE_H
#define INC_BOXTYPE_H
/// BoxType
/// Hold definitions for the 3 possible box types:
///   NONE: No box.
///   ORTHO: Orthogonal box (all angles 90.0 degrees).
///   NONORTHO: Non-orthogonal box (triclinic or truncated octahedron).
enum BoxType { NOBOX=0, ORTHO, NONORTHO };
BoxType CheckBoxType(double *,int);
BoxType SetBoxInfo(double *, double *, int);
#endif
