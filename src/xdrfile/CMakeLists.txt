#AMBER CMake buildfile for xdrfile.
set(XDRFILE_SOURCES xdrfile.c xdr_seek.c xdrfile_trr.c xdrfile_xtc.c)

#------------------------------------------------------------------------------------------
#I assume that this library should be optimized
add_compile_options(${OPT_CFLAGS})

add_library(xdrfile ${XDRFILE_SOURCES})
target_include_directories(xdrfile PUBLIC .)
install_libraries(xdrfile)