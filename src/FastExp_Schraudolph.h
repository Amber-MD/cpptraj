/*! \file FastExp_Schraudolph.h
    \brief Implementation of fast exponential by Nicol N. Schraudolph, Neural Computation 11, 853–862 (1999) 
 */
#include <math.h>
static union
{
double d;
struct
{
//#ifdef LITTLE_ENDIAN
int j, i;
//#else
//int i, j;
//#endif
} n;
} _eco;
#define EXP_A (1048576/M_LN2) /* use 1512775 for integer version */
#define EXP_C 60801 /* see text for choice of c values */
#define FASTEXPS(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)
