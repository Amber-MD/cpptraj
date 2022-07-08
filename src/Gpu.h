#ifndef INC_CPPTRAJ_GPU_H
#define INC_CPPTRAJ_GPU_H
namespace CpptrajGpu {
/// Set the major CUDA compute version number
void SetComputeVersion(int);
/// \return Max block dimensions if using 2D blocks
unsigned int MaxBlockDim_2D();
}
#endif
