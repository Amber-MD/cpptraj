#ifndef INC_XOSHIRO128PLUSPLUS_H
#define INC_XOSHIRO128PLUSPLUS_H
#ifdef __cplusplus
extern "C" {
#endif

void x128pp_seed(int);

unsigned int x128pp_next();

#ifdef __cplusplus
}
#endif
#endif
