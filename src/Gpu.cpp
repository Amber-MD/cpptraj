#include "Gpu.h"

static unsigned int cpptraj_gpu_major_ = 0;

void CpptrajGpu::SetComputeVersion(int major) {
  if (major > -1)
    cpptraj_gpu_major_ = major;
}

unsigned int CpptrajGpu::MaxBlockDim_2D() {
  if (cpptraj_gpu_major_ < 4)
    return 16;
  else
    return 32;
}
