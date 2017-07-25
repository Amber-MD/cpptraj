cmake_minimum_required(VERSION 3.1) #need 3.1 for "cmake -E env" command
project(cpptraj)

#version number
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
set(cpptraj_MAJOR_VERSION 17)
set(cpptraj_MINOR_VERSION 0)
set(cpptraj_TWEAK_VERSION 0)

set(cpptraj_VERSION "${cpptraj_MAJOR_VERSION}.${cpptraj_MINOR_VERSION}.${cpptraj_TWEAK_VERSION}")

if(NOT INSIDE_AMBER)
	# initialization and include paths
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	include(cmake/AmberBuildSystemInit.cmake)

	#now enable the languages
	enable_language(CXX)

	include(AmberBuildSystem2ndInit)
	
	include(OpenMPConfig)
	include(MPIConfig)
	include(CudaConfig)

	set(LINALG_LIBS_REQUIRED TRUE)
	set(NEEDED_3RDPARTY_TOOLS zlib fftw libbz2 readline pnetcdf blas lapack mkl arpack netcdf)
	set(REQUIRED_3RDPARTY_TOOLS netcdf)
	include(3rdPartyTools)
	
	if(fftw_DISABLED)
		enable_language(Fortran)
	endif()
	
	include(CompilerFlags)

	include(NetlibConfig)
endif()
#NetCDF and FFTW3 are found by the top-level configure script

add_subdirectory(src)

if(NOT INSIDE_AMBER)
	print_build_report()
endif()