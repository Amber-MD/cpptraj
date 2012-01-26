#ifndef INC_BOXTYPE_H
#define INC_BOXTYPE_H
/*! \file BoxType.h
    \brief Hold definitions and functions related to box information.
 */
/// High-precision value for box angle of truncated octahedron
// NOTE: Should this be in Constants.h?
#define TRUNCOCTBETA 109.4712206344906917365733534097672
/** Three possible box types:
  * - NOBOX: No box.
  * - ORTHO: Orthogonal box (all angles 90.0 degrees).
  * - NONORTHO: Non-orthogonal box (triclinic or truncated octahedron).
  */
enum BoxType { NOBOX=0, ORTHO, NONORTHO };
BoxType CheckBoxType(double *,int);
int AmberIfbox(double);
BoxType SetBoxInfo(double *, double *, int);
#endif
