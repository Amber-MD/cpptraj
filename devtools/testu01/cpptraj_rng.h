#ifndef INC_CPPTRAJ_RNG_H
#define INC_CPPTRAJ_RNG_H

#ifdef __cplusplus
extern "C" {
#endif
/// Call with type, seed
void* get_cpptraj_rng(int,int);

void destroy_cpptraj_rng(void*);

unsigned int num_cpptraj_rng(void* );
#ifdef __cplusplus
}
#endif

#endif
