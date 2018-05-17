cmake_minimum_required(VERSION 3.1) #need 3.1 for "cmake -E env" command
project(cpptraj NONE)

#version number
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
set(cpptraj_MAJOR_VERSION 4)
set(cpptraj_MINOR_VERSION 3)
set(cpptraj_TWEAK_VERSION 4)

set(cpptraj_VERSION "${cpptraj_MAJOR_VERSION}.${cpptraj_MINOR_VERSION}.${cpptraj_TWEAK_VERSION}")

if(NOT INSIDE_AMBER)
	# initialization and include paths
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	include(cmake/AmberBuildSystemInit.cmake)

	#now enable the languages
	enable_language(C CXX)

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
	
	# header installation option
	option(INSTALL_HEADERS "Copy headers to the include/cpptraj folder of the install directory.  Useful for building with pytraj." FALSE)
	
	option(BUILD_PARALLEL_COMBINATIONS "If true, then combinations of all enabled parallelizations will be built, e.g. cpptraj.OMP.MPI and cpptraj.OMP.MPI.cuda" FALSE)
else()
	set(INSTALL_HEADERS FALSE)
	set(BUILD_PARALLEL_COMBINATIONS FALSE)
endif()

# on Windows, make MSVC auto-create import libraries just like MinGW does
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

add_subdirectory(src)
add_subdirectory(test)

#--------------------------------------------------------------	
if(NOT INSIDE_AMBER)
	print_build_report()
endif()