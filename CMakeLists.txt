cmake_minimum_required(VERSION 3.3) # need 3.3 for POLICY CMP0057
project(cpptraj NONE)

#version number
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/Version.h CPPTRAJ_VERSION_FILE)
foreach(LINE ${CPPTRAJ_VERSION_FILE})
  separate_arguments(LINE)
  #message( "Line=  ${LINE}" )
  list(GET LINE 0 col0)
  if(col0 STREQUAL "#define")
    list(GET LINE 1 col1)
    if(col1 STREQUAL "CPPTRAJ_INTERNAL_VERSION")
      #message( "Col0= ${LINE}" )
      list(GET LINE 2 cpptraj_VERSION)
      break()
    endif()
  endif()
endforeach()

if(NOT INSIDE_AMBER)
	# initialization and include paths
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	include(cmake-cpptraj/BuildSystemInit.cmake)

	#now enable the languages
	enable_language(C CXX)

	include(BuildSystemSetup)

	# ARPACK support 
	option(USE_ARPACK "Whether to use the Arnoldi package for diagonalization of sparse matrices." TRUE)
	
	include(OpenMPConfig)
	include(MPIConfig)
	include(CudaConfig)

	set(LINALG_LIBS_REQUIRED TRUE)
	set(NEEDED_3RDPARTY_TOOLS zlib fftw libbz2 readline pnetcdf blas lapack mkl arpack netcdf tng_io)
        set(BUNDLED_3RDPARTY_TOOLS readline arpack tng_io)
	set(REQUIRED_3RDPARTY_TOOLS netcdf)
	include(SetupThirdParty)
	
	if((fftw_DISABLED AND mkl_DISABLED) OR arpack_ENABLED)
          # If no FFTW and no MKL we need PubFFT, so need fortran.
          # If we are building arpack we need fortran.
		enable_language(Fortran)
	endif()
	
	include(CompilerFlags)

	include(NetlibConfig)
	
	#CPack setup
	# --------------------------------------------------------------------
	set(PACKAGE_NAME "cpptraj")
	set(PACKAGE_FILENAME "cpptraj")

	#set(ICO_ICON ${CMAKE_SOURCE_DIR}/amber.ico)
	#set(ICO_UNINSTALL_ICON ${CMAKE_SOURCE_DIR}/amber-uninstall.ico)
	#set(ICNS_ICON ${CMAKE_SOURCE_DIR}/cmake-packaging/amber_logo.icns)
	#set(OSX_STARTUP_SCRIPT ${CMAKE_SOURCE_DIR}/cmake-packaging/osx-startup-script.sh)
	set(BUNDLE_IDENTIFIER org.ambermd.cpptraj)
	set(BUNDLE_SIGNATURE CPTJ)
	include(Packaging)
	
	# build options
	option(INSTALL_HEADERS "Copy headers to the include/cpptraj folder of the install directory.  Useful for building with pytraj." FALSE)
	
	option(BUILD_PARALLEL_COMBINATIONS "If true, then combinations of all enabled parallelizations will be built, e.g. cpptraj.OMP.MPI and cpptraj.OMP.MPI.cuda" FALSE)

	option(INSTALL_TESTS "Copy tests to the test/ folder of the install directory" FALSE)
else()
	set(INSTALL_HEADERS FALSE)
	set(BUILD_PARALLEL_COMBINATIONS FALSE)
endif()

# Ensure the flag for PIC is added for NVCC
if (CUDA AND CUDA_FOUND)
  if (SHARED)
    if(${COMPILER} STREQUAL GNU OR ${COMPILER} STREQUAL CLANG)
      set(CPPTRAJ_NVCC_PIC_FLAG "--compiler-options -fPIC")
    else()
      set(CPPTRAJ_NVCC_PIC_FLAG "--compiler-options -fpic")
    endif()
    message(STATUS "Set CUDA PIC flag ${CPPTRAJ_NVCC_PIC_FLAG}")
    list(APPEND CUDA_NVCC_FLAGS ${CPPTRAJ_NVCC_PIC_FLAG})
  else()
    message(STATUS "Not setting CUDA PIC flag")
  endif()
endif()

# on Windows, make MSVC auto-create import libraries just like MinGW does
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(dat)

#--------------------------------------------------------------	
if(NOT INSIDE_AMBER)
	print_build_report()
	#dump_cmake_variables()
endif()
