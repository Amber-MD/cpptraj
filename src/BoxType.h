#ifndef INC_BOXTYPE_H
#define INC_BOXTYPE_H
/// BoxType
/// Hold definitions for the 3 possible box types:
///   NOBOX: No box.
///   ORTHO: Orthogonal box (all angles 90.0 degrees).
///   NONORTHO: Non-orthogonal box (triclinic or truncated octahedron).
#define TRUNCOCTBETA 109.4712206344906917365733534097672
enum BoxType { NOBOX=0, ORTHO, NONORTHO };
BoxType CheckBoxType(double *,int);
int AmberIfbox(double);
BoxType SetBoxInfo(double *, double *, int);
#endif
