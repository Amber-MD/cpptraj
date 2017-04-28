#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
  int size = 1;
  int rank = 0;
# ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
# endif
  if (rank == 0) printf("%i\n", size);
  return 0;
}
