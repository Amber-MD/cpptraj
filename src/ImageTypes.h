#ifndef INC_IMAGETYPES_H
#define INC_IMAGETYPES_H
/*! \file ImageTypes.h
    \brief Data types and enumerations used by imaging routines. 
 */
namespace Image {
  enum Mode { BYMOL = 0, BYRES, BYATOM };
  inline const char* ModeString(Mode m) {
    if      (m == BYMOL) return "molecule";
    else if (m == BYRES) return "residue";
    else                 return "atom"; // BYATOM
  }
}
#endif
