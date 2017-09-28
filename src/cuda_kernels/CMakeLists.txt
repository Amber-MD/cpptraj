set(CPPTRAJ_CUDA_SOURCES core_kernels.cu kernel_wrappers.cu)

cuda_add_library(cpptraj_cuda_routines ${CPPTRAJ_CUDA_SOURCES})

install_libraries(cpptraj_cuda_routines)