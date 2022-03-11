#CPPTRAJ CMake buildfile for xdrfile.
set(XDRFILE_SOURCES xdrfile.c xdr_seek.c xdrfile_trr.c xdrfile_xtc.c)

#------------------------------------------------------------------------------------------
# Optimize library
add_compile_options(${OPT_CFLAGS})

# Always build bundled xdrfile as static
add_library(xdrfile STATIC ${XDRFILE_SOURCES})
make_pic_if_needed( xdrfile )

target_include_directories(xdrfile PUBLIC .)
install_libraries(xdrfile)
