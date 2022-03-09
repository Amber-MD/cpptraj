#CPPTRAJ CMake buildfile for arpack.
set(ARPACK_SOURCES ivout.f dvout.f dgetv0.f dlaqrb.f dnaitr.f dnapps.f 
dnaup2.f dnaupd.f dnconv.f dneigh.f dmout.f 
dneupd.f dngets.f dsaitr.f dsapps.f dsaup2.f 
dsaupd.f dsconv.f dseigt.f dsesrt.f dseupd.f 
dsgets.f dsortc.f dsortr.f dstatn.f dstats.f dstqrb.f arsecond.f)

#------------------------------------------------------------------------------------------
#I assume that this library should be optimized
add_compile_options(${OPT_FFLAGS})

# Add shared flag if needed
#foreach(ARPACK_SOURCE_FILE ${ARPACK_SOURCES})
#  set_property(SOURCE ${ARPACK_SOURCE_FILE} PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})
#endforeach()

#if("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
#	# xdrfile is not set up to build as a Windows shared library.  Someone would need to add dll exports to all of the symbols (bluh)
#	add_library(xdrfile STATIC ${XDRFILE_SOURCES})
#else()
	add_library(arpack STATIC ${ARPACK_SOURCES})
        make_pic_if_needed( arpack )
#endif()

target_include_directories(arpack PUBLIC .)
install_libraries(arpack)
