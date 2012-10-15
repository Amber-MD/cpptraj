#ifndef INC_THERMO_H
#define INC_THERMO_H
#include "CpptrajFile.h"
/*! \file Thermo.h
    \brief Calculate thermochemistry. 
 */
void thermo( CpptrajFile&, int, int, int, const double*, const double*, 
             const double*, double, double);
#endif
