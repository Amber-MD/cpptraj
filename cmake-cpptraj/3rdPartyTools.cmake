#  This file configures which 3rd party tools are built in and which are used from the system.
#  NOTE: must be included after WhichTools.cmake

message(STATUS "Checking whether to use built-in libraries...")


#List of 3rd party tools.
set(3RDPARTY_TOOLS
blas
lapack
arpack 
readline
netcdf
pnetcdf
fftw
zlib
libbz2
libm
mkl
tng_io)

set(3RDPARTY_TOOL_USES 
"for fundamental linear algebra calculations"                                     
"for fundamental linear algebra calculations"                                     
"for fundamental linear algebra calculations"                                     
"enables an interactive terminal"
"for creating NetCDF trajectory data files"                                              
"for creating NetCDF trajectory files in parallel"                                  
"used to do Fourier transforms very quickly"                                      
"for various compression and decompression tasks"                                 
"for bzip2 compression in cpptraj"  
"for fundamental math routines if they are not contained in the C library"                                              
"alternate implementation of lapack and blas that is tuned for speed"             
"enables GROMACS tng trajectory input in cpptraj")                                                  

# Logic to disable tools
set(3RDPARTY_SUBDIRS "")

#sets a tool to external, internal, or disabled
#STATUS=EXTERNAL, INTERNAL, or DISABLED
macro(set_3rdparty TOOL STATUS)
	set(${TOOL}_INTERNAL FALSE)
	set(${TOOL}_EXTERNAL FALSE)
	set(${TOOL}_DISABLED FALSE)
	set(${TOOL}_ENABLED TRUE)
	
	
	if(${STATUS} STREQUAL EXTERNAL)
		set(${TOOL}_EXTERNAL TRUE)
	elseif(${STATUS} STREQUAL INTERNAL)
		
		#the only way to get this message would be to use FORCE_INTERNAL_LIBS incorrectly, unless someone messed up somewhere
		if("${BUNDLED_3RDPARTY_TOOLS}" MATCHES ${TOOL})
				set(${TOOL}_INTERNAL TRUE)
		else()
			if(INSIDE_AMBER)
				# getting here means there's been a programming error
				message(FATAL_ERROR "3rd party program ${TOOL} is not bundled and cannot be built inside Amber.")
			elseif("${REQUIRED_3RDPARTY_TOOLS}" MATCHES ${TOOL})
				# kind of a kludge - even when we're in a submodule, things will still get set to internal, so we just treat internal equal to disabled
				message(FATAL_ERROR "3rd party program ${TOOL} is required, but was not found.")
			else()
				# we're in a submodule, and it's not required, so it's OK that the tool is not bundled
				set(${TOOL}_DISABLED TRUE)
				set(${TOOL}_ENABLED FALSE)
				
			endif()
		endif()
	
	else()
		list_contains(TOOL_REQUIRED ${TOOL} ${REQUIRED_3RDPARTY_TOOLS})
				
		if(TOOL_REQUIRED)
			message(FATAL_ERROR "3rd party program ${TOOL} is required to build ${PROJECT_NAME}, but it is disabled (likely because it was not found).")
		endif()
		
		set(${TOOL}_DISABLED TRUE)
		set(${TOOL}_ENABLED FALSE)
		
	endif()	
endmacro(set_3rdparty)

#------------------------------------------------------------------------------
#  OS threading library (not really a 3rd party tool)
#------------------------------------------------------------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE) #Yeah, we're biased.
find_package(Threads)

# first, figure out which tools we need
# -------------------------------------------------------------------------------------------------------------------------------

# if NEEDED_3RDPARTY_TOOLS is not passed in, assume that all of them are needed
if(NOT DEFINED NEEDED_3RDPARTY_TOOLS)
	set(NEEDED_3RDPARTY_TOOLS "${3RDPARTY_TOOLS}")
endif()

if(NOT DEFINED BUNDLED_3RDPARTY_TOOLS)
	set(BUNDLED_3RDPARTY_TOOLS "")
endif()

# suspicious 3rd party tools: tools where there are often problems with the system install
if(NOT DEFINED SUSPICIOUS_3RDPARTY_TOOLS)
	set(SUSPICIOUS_3RDPARTY_TOOLS "mkl,boost")
endif()

option(TRUST_SYSTEM_LIBS "If true, Amber will use all found system libraries, even if they are considered problematic by the Amber developers" FALSE)

foreach(TOOL ${3RDPARTY_TOOLS})
	list(FIND NEEDED_3RDPARTY_TOOLS ${TOOL} TOOL_INDEX)
	
	test(NEED_${TOOL} NOT "${TOOL_INDEX}" EQUAL -1)
endforeach()

if(("${NEEDED_3RDPARTY_TOOLS}" MATCHES "mkl" OR "${NEEDED_3RDPARTY_TOOLS}" MATCHES "blas" OR "${NEEDED_3RDPARTY_TOOLS}" MATCHES "lapack")
	AND NOT ("${NEEDED_3RDPARTY_TOOLS}" MATCHES "mkl" AND "${NEEDED_3RDPARTY_TOOLS}" MATCHES "blas" AND "${NEEDED_3RDPARTY_TOOLS}" MATCHES "lapack"))
	message(FATAL_ERROR "If any of mkl, blas, and lapack are put into NEEDED_3RDPARTY_TOOLS, them you must put all of them in since mkl replaces blas and lapack")
endif()

# 1st pass checking
# -------------------------------------------------------------------------------------------------------------------------------

# Force all try_link_library() calls to build position independent executables if Amber is built as shared.
# This ensures that all dependent libraries are built as position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ${SHARED})

#------------------------------------------------------------------------------
#  Readline
#------------------------------------------------------------------------------

if(NEED_readline)
	find_package(Readline)
	
	if(READLINE_FOUND)
		set_3rdparty(readline EXTERNAL)
	elseif(NOT TARGET_WINDOWS)
		set_3rdparty(readline INTERNAL)
	else()
		set_3rdparty(readline DISABLED)
	endif()
endif()

#------------------------------------------------------------------------------
#  MKL (near the top because it contains lapack, and blas)
#------------------------------------------------------------------------------

if(NEED_mkl)

	set(MKL_NEEDINCLUDES FALSE)
	set(MKL_NEEDEXTRA FALSE)
	
	# Static MKL is not the default at this time.
	# <long_explanation>
	# MKL has a fftw3 compatibility interface.  Wierdly enough, this interface is spread out between several different libraries: the main interface library, the 
	# cdft library, and the actual fftw3 interface library (which is distributed as source code, not a binary).
	# So, even though we don't use the fftw3 interface, there are symbols in the main MKL libraries which conflict with the symbols from fftw3.
	# Oddly, on many platforms, the linker handles this fine.  However, in at least one case (the SDSC supercomputer Comet, running a derivative of CentOS),
	# ld balks at this multiple definition, and refuses to link programs which use MKL and fftw, but ONLY when BOTH of them are built as static libraries.
	# Why this is, I'm not sure.  I do know that it's better to build fftw3 as static and use mkl as shared (because mkl is a system library)
	# then the other way around, so that's what I do
	# </long_explanation>
	option(MKL_STATIC "Whether to prefer MKL's static libraries" FALSE)

	option(MKL_MULTI_THREADED "Whether to link MKL in OpenMP mode to parallelize singlethreaded calls to MKL functions" TRUE)
	
	set(MKL_STATIC FALSE)
	find_package(MKL)
	
	if(MKL_FOUND)
		set_3rdparty(mkl EXTERNAL)
	else()
		set_3rdparty(mkl DISABLED)
	endif()
endif()


#------------------------------------------------------------------------------
#  FFTW
#------------------------------------------------------------------------------
if(NEED_fftw)

	if(DEFINED USE_FFT AND NOT USE_FFT)
		set_3rdparty(fftw DISABLED)
	else()	
          find_package(FFTW)	
		#if(MPI)
		#	find_package(FFTW COMPONENTS MPI Fortran)
		#else()
		#	find_package(FFTW COMPONENTS Fortran)
		#endif()
	
		if(FFTW_FOUND)
			set_3rdparty(fftw EXTERNAL)
		else()
			set_3rdparty(fftw INTERNAL)
		endif()
	endif()
endif()

#------------------------------------------------------------------------------
#  NetCDF
#------------------------------------------------------------------------------

if(NEED_netcdf OR NEED_netcdf-fortran)
	#find_package(NetCDF OPTIONAL_COMPONENTS F77 F90)
	find_package(NetCDF)
endif()

if(NEED_netcdf)
	if(NetCDF_FOUND)
		set_3rdparty(netcdf EXTERNAL)
	else()
		set_3rdparty(netcdf INTERNAL)
	endif()
endif()

if(NEED_netcdf-fortran)
	if(NetCDF_F90_FOUND)
		set_3rdparty(netcdf-fortran EXTERNAL)
	else()
		set_3rdparty(netcdf-fortran INTERNAL)
	endif()
endif()

#------------------------------------------------------------------------------
#  Netlib libraries
#------------------------------------------------------------------------------

if(NEED_blas) # because of the earlier check, we can be sure that NEED_blas == NEED_lapack

	# this calls FindBLAS
	find_package(LAPACKFixed)
	
	if(BLAS_FOUND)
		set_3rdparty(blas EXTERNAL)
	else()
		set_3rdparty(blas INTERNAL)
	endif()
	
	if(LAPACK_FOUND)
		set_3rdparty(lapack EXTERNAL)
	else()
		set_3rdparty(lapack INTERNAL)
	endif()
endif()

if(NEED_arpack)
	#  ARPACK
	find_package(ARPACK)
	
	if(ARPACK_FOUND)
		set_3rdparty(arpack EXTERNAL)
	else()
		set_3rdparty(arpack INTERNAL)
	endif()
endif()

# --------------------------------------------------------------------
#  Parallel NetCDF
# --------------------------------------------------------------------

if(NEED_pnetcdf)
	find_package(PnetCDF)
	
	if(PnetCDF_FOUND)
		set_3rdparty(pnetcdf EXTERNAL)
	else()
		set_3rdparty(pnetcdf INTERNAL)
	endif()
endif()

#------------------------------------------------------------------------------
#  zlib, for cpptraj and netcdf
#------------------------------------------------------------------------------

if(NEED_zlib)
	find_package(ZLIB)
	
	if(ZLIB_FOUND)
		set_3rdparty(zlib EXTERNAL)
	else()
		set_3rdparty(zlib DISABLED)  # will always error
	endif()
endif()

#------------------------------------------------------------------------------
#  bzip2
#------------------------------------------------------------------------------
if(NEED_libbz2)
	find_package(BZip2)
	
	
	if(BZIP2_FOUND)
		set_3rdparty(libbz2 EXTERNAL)
	else()
		set_3rdparty(libbz2 DISABLED)
	endif()
endif()

#------------------------------------------------------------------------------
#  Math library
#------------------------------------------------------------------------------ 

if(NEED_libm)
	# figure out if we need a math library
	find_package(CMath)

	if(CMath_FOUND)
		set_3rdparty(libm EXTERNAL)
	else()
		set_3rdparty(libm DISABLED)
	endif()
endif()

#------------------------------------------------------------------------------
#  tng_io
#------------------------------------------------------------------------------

if(NEED_tng_io)
	find_package(tng_io CONFIG)

	if(NOT tng_io_FOUND)
		message(STATUS "Could not find tng_io.  To locate it, add its install dir to the prefix path.")
	endif()

	if(tng_io_FOUND)
		set_3rdparty(tng_io EXTERNAL)
	elseif(zlib_ENABLED) # tng io requires zlib
		set_3rdparty(tng_io INTERNAL)
	else()
		set_3rdparty(tng_io DISABLED)
	endif() 
endif()

# Apply overrides
# -------------------------------------------------------------------------------------------------------------------------------------------------------

# we can now reset this back to the default behavior -- targets will be made PIC as needed in the individual CMake scripts
unset(CMAKE_POSITION_INDEPENDENT_CODE)

set(FORCE_EXTERNAL_LIBS "" CACHE STRING "3rd party libraries to force using the system version of. Accepts a semicolon-seperated list of library names from the 3rd Party Libraries section of the build report.")
set(FORCE_INTERNAL_LIBS "" CACHE STRING "3rd party libraries to force to build inside Amber. Accepts a semicolon-seperated list of library names from the 3rd Party Libraries section of the build report.")
set(FORCE_DISABLE_LIBS "" CACHE STRING "3rd party libraries to force Amber to not use at all. Accepts a semicolon-seperated list of library names from the 3rd Party Libraries section of the build report.")

# look for and handle suspicious tools
if(NOT TRUST_SYSTEM_LIBS)
	set(FOUND_SUSPICIOUS_TOOLS FALSE)
	set(SUSPICIOUS_TOOLS_MESSAGE "")

	foreach(TOOL ${SUSPICIOUS_3RDPARTY_TOOLS})
		if(${TOOL}_EXTERNAL)
			set(FOUND_SUSPICIOUS_TOOLS TRUE)
			set(SUSPICIOUS_TOOLS_MESSAGE "${SUSPICIOUS_TOOLS_MESSAGE} ${TOOL}")

			if("${TOOL}" IN_LIST BUNDLED_3RDPARTY_TOOLS)
				set_3rdparty(${TOOL} INTERNAL)
			else()
				set_3rdparty(${TOOL} DISABLED)
			endif()
		endif()
	endforeach()

	if(FOUND_SUSPICIOUS_TOOLS)
		message(STATUS "The following libraries were found on your system, but are not being used:${SUSPICIOUS_TOOLS_MESSAGE}")
		message(STATUS "This is because the Amber devs have frequently seen broken installs of these libraries causing trouble.")
		message(STATUS "To tell ${PROJECT_NAME} to link to these libraries, you can select individual ones using -DFORCE_EXTERNAL_LIBS=<names>.")
		message(STATUS "Alternatively, you can request that ${PROJECT_NAME} link to all libraries it finds using -DTRUST_SYSTEM_LIBS=TRUE.")
	endif()
endif()

foreach(TOOL ${FORCE_EXTERNAL_LIBS})
	colormsg(YELLOW "Forcing ${TOOL} to be sourced externally")

	list_contains(VALID_TOOL ${TOOL} ${3RDPARTY_TOOLS})
	
	if(NOT VALID_TOOL)
		message(FATAL_ERROR "${TOOL} is not a valid 3rd party library name.")
	endif()
	
	set_3rdparty(${TOOL} EXTERNAL)
endforeach()

if(INSIDE_AMBER)
	foreach(TOOL ${FORCE_INTERNAL_LIBS})
		colormsg(GREEN "Forcing ${TOOL} to be built internally")
	
		list_contains(VALID_TOOL ${TOOL} ${3RDPARTY_TOOLS})
		
		if(NOT VALID_TOOL)
			message(FATAL_ERROR "${TOOL} is not a valid 3rd party library name.")
		endif()
		
		set_3rdparty(${TOOL} INTERNAL)
	endforeach()
endif()

foreach(TOOL ${FORCE_DISABLE_LIBS})
	colormsg(HIRED "Forcing ${TOOL} to be disabled")

	list_contains(VALID_TOOL ${TOOL} ${3RDPARTY_TOOLS})
	
	if(NOT VALID_TOOL)
		message(FATAL_ERROR "${TOOL} is not a valid 3rd party library name.")
	endif()
	
	set_3rdparty(${TOOL} DISABLED)
endforeach()

# force all unneeded tools to be disabled
foreach(TOOL ${3RDPARTY_TOOLS})
	list(FIND NEEDED_3RDPARTY_TOOLS ${TOOL} TOOL_INDEX)
	
	if(${TOOL_INDEX} EQUAL -1)
		set_3rdparty(${TOOL} DISABLED)
	endif()
	
endforeach()

# check math library configuration
if(LINALG_LIBS_REQUIRED AND NOT (mkl_ENABLED OR (blas_ENABLED AND lapack_ENABLED)))
	message(FATAL_ERROR "You must enable a linear algebra library -- either blas and lapack, or mkl")
endif()

if(mkl_ENABLED AND (blas_ENABLED AND lapack_ENABLED))
	# prefer MKL to BLAS
	set_3rdparty(blas DISABLED)
	set_3rdparty(lapack DISABLED)
endif()

# Now that we know which libraries we need, set them up properly.
# -------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#  Readline
#------------------------------------------------------------------------------

if(readline_EXTERNAL)

	# rename target
	import_libraries(readline 
		LIBRARIES readline::readline 
		INCLUDES ${READLINE_INCLUDE_DIR}/readline)	
elseif(readline_INTERNAL)

	list(APPEND 3RDPARTY_SUBDIRS cpptraj/src/readline)

endif()

#------------------------------------------------------------------------------
#  MKL
#------------------------------------------------------------------------------


if(mkl_ENABLED)
	if(NOT MKL_FOUND)
		message(FATAL_ERROR "You enabled MKL, but it was not found.")
	endif()
	
	if(MIXING_COMPILERS AND OPENMP)
		message(WARNING "You are using different compilers from different vendors together.  This may cause link errors with MKL and OpenMP.  There is no way around this.")
	endif()
	
	if(mkl_ENABLED AND (blas_ENABLED OR lapack_ENABLED))
		message(FATAL_ERROR "MKL replaces blas and lapack!  They can't be enabled when MKL is in use!")
	endif()
	
	# add to library tracker
	import_libraries(mkl LIBRARIES ${MKL_FORTRAN_LIBRARIES} INCLUDES ${MKL_INCLUDE_DIRS})
endif()

#------------------------------------------------------------------------------
#  FFTW
#------------------------------------------------------------------------------

if(fftw_EXTERNAL)
	
	if(NOT FFTW_FOUND)
		message(FATAL_ERROR "Could not find FFTW, but it was set so be sourced externally!")
	endif()

	# rename targets
	import_libraries(fftw LIBRARIES fftw::fftwl)

	if(MPI)
		import_libraries(fftw_mpi LIBRARIES fftw::mpi_l)
	endif()	
	
elseif(fftw_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS fftw-3.3)
endif()


#------------------------------------------------------------------------------
#  NetCDF
#------------------------------------------------------------------------------

if(netcdf_EXTERNAL)
	
	if(NOT NetCDF_FOUND)
		message(FATAL_ERROR "netcdf was set to be sourced externally, but it was not found!")
	endif()
	
	# Import the system netcdf as a library
	import_library(netcdf ${NetCDF_LIBRARIES_C} ${NetCDF_INCLUDES})
	
elseif(netcdf_INTERNAL)

	if(${COMPILER} STREQUAL CRAY)
			message(FATAL_ERROR "Bundled NetCDF cannot be used with cray compilers.  Please reconfigure with -DFORCE_EXTERNAL_LIBS=netcdf. \
		 On cray systems you can usually load the system NetCDF with 'module load cray-netcdf' or 'module load netcdf'.")
	endif()
		
	list(APPEND 3RDPARTY_SUBDIRS netcdf-4.6.1)
endif()

if(netcdf-fortran_EXTERNAL)

	if(NOT NetCDF_F90_FOUND)
		message(FATAL_ERROR "netcdf-fortran was set to be sourced externally, but it was not found!")
	endif()

	if(netcdf_INTERNAL)
		message(FATAL_ERROR "Cannot use internal netcdf with external netcdf-fortran!")
	endif()
	
	# Import the system netcdf as a library
	import_library(netcdff ${NetCDF_LIBRARIES_F90} ${NetCDF_INCLUDES_F90})
	set_property(TARGET netcdff PROPERTY INTERFACE_LINK_LIBRARIES netcdf)
	
	# This is really for symmetry with the other MOD_DIRs more than anything.
	set(NETCDF_FORTRAN_MOD_DIR ${NetCDF_INCLUDES_F90})
		
elseif(netcdf-fortran_INTERNAL)
	#TODO on Cray systems a static netcdf may be required

	if(${COMPILER} STREQUAL CRAY)
			message(FATAL_ERROR "Bundled NetCDF cannot be used with cray compilers. \
 On cray systems you can usually load the system NetCDF with 'module load cray-netcdf' or 'module load netcdf'.")
	endif()

	list(APPEND 3RDPARTY_SUBDIRS netcdf-fortran-4.4.4)
endif()

#------------------------------------------------------------------------------
#  Netlib libraries
#------------------------------------------------------------------------------

# BLAS
if(blas_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS blas)
elseif(blas_EXTERNAL)
	import_libraries(blas LIBRARIES "${BLAS_LIBRARIES}")
endif()

#  LAPACK
if(lapack_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS lapack)
elseif(lapack_EXTERNAL)
	import_libraries(lapack LIBRARIES "${LAPACK_LIBRARIES}")
endif()


#  ARPACK
if(arpack_EXTERNAL)
	if(NOT EXISTS "${ARPACK_LIBRARY}")
		message(FATAL_ERROR "arpack was set to be sourced externally, but it was not found!")
	endif()
	
	import_library(arpack ${ARPACK_LIBRARY})

	if(NOT ARPACK_HAS_ARSECOND)
		message(STATUS "System arpack is missing the arsecond_ function.  That function will be built inside amber")
	endif()

elseif(arpack_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS arpack)
endif()
	
# --------------------------------------------------------------------
#  Parallel NetCDF
# --------------------------------------------------------------------

if(pnetcdf_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS pnetcdf)
elseif(pnetcdf_EXTERNAL)
	if(NOT PnetCDF_FOUND)
		message(FATAL_ERROR "You requested to use an external pnetcdf, but no installation was found.")
	endif()
	
	# rename target
	import_libraries(pnetcdf LIBRARIES pnetcdf::pnetcdf)
	
endif()

#------------------------------------------------------------------------------
# libbz2
#------------------------------------------------------------------------------

if(libbz2_EXTERNAL)
	using_library_targets(BZip2::BZip2)
endif()

# --------------------------------------------------------------------
#  tng_io
# --------------------------------------------------------------------

if(tng_io_INTERNAL)
	list(APPEND 3RDPARTY_SUBDIRS cpptraj/src/tng)
elseif(tng_io_EXTERNAL)
	if(NOT tng_io_FOUND)
		message(FATAL_ERROR "You requested to use an external libtng_io, but no installation was found.")
	endif()
	
	# rename target
	import_libraries(tng_io LIBRARIES tng_io::tng_io)
endif() 
