#CPPTRAJ CMake buildfile for arpack.
set(ARPACK_SOURCES ivout.o dvout.o dgetv0.o dlaqrb.o dnaitr.o dnapps.o 
dnaup2.o dnaupd.o dnconv.o dneigh.o dmout.o 
dneupd.o dngets.o dsaitr.o dsapps.o dsaup2.o 
dsaupd.o dsconv.o dseigt.o dsesrt.o dseupd.o 
dsgets.o dsortc.o dsortr.o dstatn.o dstats.o dstqrb.o arsecond.o)

#------------------------------------------------------------------------------------------
#I assume that this library should be optimized
add_compile_options(${OPT_FFLAGS})

#if("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
#	# xdrfile is not set up to build as a Windows shared library.  Someone would need to add dll exports to all of the symbols (bluh)
#	add_library(xdrfile STATIC ${XDRFILE_SOURCES})
#else()
	add_library(arpack ${ARPACK_SOURCES})
#endif()

target_include_directories(arpack PUBLIC .)
install_libraries(arpack)
