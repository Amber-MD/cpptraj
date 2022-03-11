#CPPTRAJ CMake buildfile for arpack.
set(ARPACK_SOURCES ivout.f dvout.f dgetv0.f dlaqrb.f dnaitr.f dnapps.f 
dnaup2.f dnaupd.f dnconv.f dneigh.f dmout.f 
dneupd.f dngets.f dsaitr.f dsapps.f dsaup2.f 
dsaupd.f dsconv.f dseigt.f dsesrt.f dseupd.f 
dsgets.f dsortc.f dsortr.f dstatn.f dstats.f dstqrb.f arsecond.f)

#------------------------------------------------------------------------------------------
# Optimize library
add_compile_options(${OPT_FFLAGS})

# Always build the bundled arpack as static.
add_library(arpack STATIC ${ARPACK_SOURCES})
make_pic_if_needed( arpack )

target_include_directories(arpack PUBLIC .)
install_libraries(arpack)
