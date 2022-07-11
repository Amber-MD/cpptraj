set(CPPTRAJ_CUDA_SOURCES core_kernels.cu kernel_wrappers.cu GistCudaCalc.cu GistCudaSetup.cu kernel_rdf.cu)

cuda_add_library(cpptraj_cuda_routines STATIC ${CPPTRAJ_CUDA_SOURCES})
#make_pic_if_needed does not appear to work for NVCC. Instead, added
# logic to $CPPTRAJHOME/CMakeLists.txt
#make_pic_if_needed(cpptraj_cuda_routines)

install_libraries(cpptraj_cuda_routines)
