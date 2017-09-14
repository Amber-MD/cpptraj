set(CPPTRAJ_CUDA_SOURCES core_kernels.cu kernel_wrappers.cu)

cuda_add_library(cpptraj_cuda STATIC ${CPPTRAJ_CUDA_SOURCES})