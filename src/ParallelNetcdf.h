#ifndef INC_PARALLELNETCDF_H
#define INC_PARALLELNETCDF_H
/** In C++, MPI is a reserved namespace, but since Amber uses it to define
  * parallel build we need to work with it. Change MPI to CPPTRAJ_MPI for
  * this header file only.
  */
#ifdef MPI
# undef MPI
# define CPPTRAJ_MPI
#endif
#ifdef CPPTRAJ_MPI
# ifdef HAS_PNETCDF
#   include <pnetcdf.h>
int checkPNCerr(int);
# endif /* HAS_PNETCDF */
# undef CPPTRAJ_MPI
# define MPI
#endif /* CPPTRAJ_MPI */
#endif
