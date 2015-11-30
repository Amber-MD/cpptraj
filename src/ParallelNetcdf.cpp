#include "ParallelNetcdf.h"
#include "Parallel.h"
#include "CpptrajStdio.h"

int checkPNCerr(int err) {
  if (err != NC_NOERR) {
    rprinterr("PNetCDF Error: %s\n", ncmpi_strerror(err));
    Parallel::Abort( err );
    return 1;
  }
  return 0;
}
