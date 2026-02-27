#ifndef INC_STRUCTURE_STRUCTUREENUM_H
#define INC_STRUCTURE_STRUCTUREENUM_H
namespace Cpptraj {
/** @brief The namespace for all structure related classes.
  *
  * This namespace holds all classes related to building systems.
  */
namespace Structure {

enum ChiralType { CHIRALITY_ERR = 0, IS_S, IS_R, IS_UNKNOWN_CHIRALITY };
/// \return String corresponding to ChiralType
const char* chiralStr(ChiralType);

/// Residue terminal type
enum TerminalType { BEG_TERMINAL = 0, NON_TERMINAL, END_TERMINAL };
/// \return String corresponding to TerminalType
const char* terminalStr(TerminalType);

}
}
#endif
