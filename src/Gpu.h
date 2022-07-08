#ifndef INC_CPPTRAJ_GPU_H
#define INC_CPPTRAJ_GPU_H
namespace Gpu {
#ifdef CUDA
void SetComputeVersion(int);
#endif
}
#endif
