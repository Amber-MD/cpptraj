#ifndef INC_MOLSURF_H
#define INC_MOLSURF_H
#ifdef __cplusplus
extern "C" {
#endif

#define REAL_T	double
#include "Name.h" // NAME

REAL_T molsurf(REAL_T, REAL_T *, int, NAME *, NAME *, int *, REAL_T *, REAL_T *);

#ifdef __cplusplus
}
#endif
#endif
