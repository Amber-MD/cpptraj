#ifdef CUDA
#include "Gpu.h"

static unsigned int cpptraj_gpu_major_ = 0;

void Gpu::SetComputeVersion(int major) {
  if (major > -1)
    cpptraj_gpu_major_ = major;
}
#endif /* CUDA */
