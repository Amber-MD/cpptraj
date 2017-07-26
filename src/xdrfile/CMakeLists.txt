#AMBER CMake buildfile for xdrfile.
set(XDRFILE_SOURCES xdrfile.c xdr_seek.c xdrfile_trr.c xdrfile_xtc.c)

#------------------------------------------------------------------------------------------
#I assume that this library should be optimized
add_compile_options(${OPT_CFLAGS})

if("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
	# xdrfile is not set up to build as a Windows shared library.  Someone would need to add dll exports to all of the symbols (bluh)
	add_library(xdrfile STATIC ${XDRFILE_SOURCES})
else()
	add_library(xdrfile ${XDRFILE_SOURCES})
endif()

target_include_directories(xdrfile PUBLIC .)
install_libraries(xdrfile)