#include "StructureEnum.h"

/** \return string corresponding to ChiralType */
const char* Cpptraj::Structure::chiralStr(ChiralType ct) {
  switch (ct) {
    case IS_S : return "S";
    case IS_R : return "R";
    case IS_UNKNOWN_CHIRALITY : return "Unknown";
  }
  return 0;
}

/** return string corresponding to TerminalType */
const char* Cpptraj::Structure::terminalStr(TerminalType tt) {
  switch (tt) {
    case BEG_TERMINAL : return "Begin";
    case NON_TERMINAL : return "Non";
    case END_TERMINAL : return "End";
  }
  return 0;
}
