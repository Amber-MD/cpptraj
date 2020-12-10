#ifndef INC_FASTEXP_SCHRAUDOLPH_H
#define INC_FASTEXP_SCHRAUDOLPH_H
/*! \file FastExp_Schraudolph.h
    \brief Implementations of fast exponential by Nicol N. Schraudolph, Neural Computation 11, 853–862 (1999) 
 */
#include <math.h> // M_LN2
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

/// Approximate exp adapted from Schraudolph, 1999 - double precision floating point version.
/** \sa https://bduvenhage.me/performance/machine_learning/2019/06/04/fast-exp.html
  */
ALWAYS_INLINE double fast_exps_64(const double x) noexcept {
    // Based on Schraudolph 1999, A Fast, Compact Approximation of the Exponential Function.
    // - Adapted to use 64-bit integer; reduces staircase effect.
    // - Valid for x in approx range (-700, 700).
    union{double d_; int64_t i_;} uid; //This could be moved to the thread scope.
    //BBBD(sizeof(uid)!=8)
    uid.i_ = int64_t(double((int64_t(1) << 52) / log(2.0)) * x + double((int64_t(1) << 52) * 1023 - 0)); //c=0 for 1.0 at zero.
    return uid.d_;
}
#endif
