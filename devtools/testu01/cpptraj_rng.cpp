#include "cpptraj_rng.h"
#include "../../src/Random.h"

void* get_cpptraj_rng(int itype, int iseed) {
  Random_Number::SetDefaultRng((Random_Number::RngType)itype);
  Random_Number* rng = new Random_Number();
  rng->rn_set( iseed );

  return ( reinterpret_cast<void*>( rng ) );
}

void destroy_cpptraj_rng(void* rng) {
  delete( reinterpret_cast<Random_Number*>( rng ) );
}

unsigned int num_cpptraj_rng(void* rngIn) {
  Random_Number* rng = reinterpret_cast<Random_Number*>( rngIn );
  return rng->rn_num();
}
