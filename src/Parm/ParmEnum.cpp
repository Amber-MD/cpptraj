#include "ParmEnum.h"

const char* Cpptraj::Parm::WaterModelStr(WaterModelType wm) {
  switch (wm) {
    case TIP3P : return "TIP3P"; 
    case TIP4PEW : return "TIP4Pew"; 
    case SPCE : return "SPCe"; 
    case OPC3 : return "OPC3"; 
    case OPC : return "OPC"; 
    case FB3 : return "FB3"; 
    case FB4 : return "FB4"; 
    case UNKNOWN_WATER_MODEL: return "Unknown water model";
  }
  return 0;
}
