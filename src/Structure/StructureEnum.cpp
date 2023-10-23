#include "StructureEnum.h"

/** \return string corresponding to ChiralType */
const char* Cpptraj::Structure::chiralStr(ChiralType ct) {
  switch (ct) {
    case CHIRALITY_ERR : return "Error";
    case IS_S : return "S";
    case IS_R : return "R";
    case IS_UNKNOWN_CHIRALITY : return "Unknown";
  }
  return 0;
}
