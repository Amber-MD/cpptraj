#ifndef INC_PARM_PARMENUM_H
#define INC_PARM_PARMENUM_H
namespace Cpptraj {
/** @brief Namespace containing classes related to parameters.
  *
  * This namespace contains classes/functions related to assigning, updating, and
  * manipulating paramters.
  */
namespace Parm {
  enum RetType { ADDED = 0, SAME, UPDATED, ERR };

  enum WaterModelType { TIP3P = 0, TIP4PEW, SPCE, OPC3, OPC, FB3, FB4, UNKNOWN_WATER_MODEL };
  const char* WaterModelStr(WaterModelType);
}
}
#endif
