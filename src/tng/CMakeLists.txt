# Amber buildfile for cpptraj's internal libtng_io.

set(TNG_IO_SOURCES
	md5.c 
	tng_io.c 
	tng_compress.c 
	coder.c 
	warnmalloc.c 
	fixpoint.c 
	xtc3.c 
	bwlzh.c 
	xtc2.c 
	widemuldiv.c 
	huffmem.c 
	vals16.c 
	bwt.c 
	mtf.c 
	rle.c 
	lz77.c 
	dict.c 
	huffman.c 
	merge_sort.c)

add_compile_options(${OPT_CFLAGS})

# add library
add_library(tng_io STATIC ${TNG_IO_SOURCES})
make_pic_if_needed(tng_io)
target_include_directories(tng_io PUBLIC .)

# checks based on tng_io's cmake script
check_all_includes(C inttypes.h)
if (TNG_HAVE_INTTYPES_H)
    target_compile_definitions(tng_io PUBLIC USE_STD_INTTYPES_H)
endif()

include(TestBigEndian)
test_big_endian(BIGENDIAN)
if(BIGENDIAN)
  set_property(SOURCE md5.c APPEND PROPERTY COMPILE_DEFINITIONS TNG_INTEGER_BIG_ENDIAN)
endif()


