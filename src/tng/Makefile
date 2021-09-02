# Makefile for TNG bundled with cpptraj
include ../../external.config.h

CFLAGS += -I . $(ZLIB_INCL)

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libtng_io.a

# Source files
SOURCES=md5.c tng_io.c tng_compress.c coder.c warnmalloc.c fixpoint.c xtc3.c bwlzh.c xtc2.c widemuldiv.c huffmem.c vals16.c bwt.c mtf.c rle.c lz77.c dict.c huffman.c merge_sort.c

# Objects
OBJECTS=$(SOURCES:.c=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

test: $(TARGET) tng_io_testing.o tng_example.tng
	$(CC) -o a.out tng_io_testing.c $(TARGET) -lm $(ZLIB_FLAG)
	./a.out

clean:
	$(DEL_FILE) *.o $(TARGET) tng_io_testing.o a.out tng_example_out.tng tng_test.tng

uninstall: clean

# Dependencies
md5.o : md5.c tng/md5.h

tng_io.o : tng_io.c tng/tng_io.h tng/tng_io_fwd.h compression/tng_compress.h tng/version.h

tng_compress.o : tng_compress.c compression/tng_compress.h compression/coder.h compression/fixpoint.h

coder.o : coder.c compression/coder.h compression/tng_compress.h compression/bwlzh.h compression/warnmalloc.h

warnmalloc.o : warnmalloc.c compression/warnmalloc.h compression/tng_compress.h

fixpoint.o : fixpoint.c compression/fixpoint.h 

xtc3.o : xtc3.c compression/warnmalloc.h compression/bwlzh.h compression/widemuldiv.h

bwlzh.o : bwlzh.c compression/bwlzh.h compression/warnmalloc.h compression/tng_compress.h compression/dict.h compression/vals16.h compression/rle.h compression/mtf.h compression/bwt.h compression/lz77.h

xtc2.o : xtc2.c compression/coder.h compression/widemuldiv.h compression/warnmalloc.h

widemuldiv.o : widemuldiv.c compression/widemuldiv.h compression/tng_compress.h compression/my64bit.h

huffmem.o : huffmem.c compression/huffman.h compression/warnmalloc.h compression/tng_compress.h compression/bwlzh.h compression/dict.h compression/rle.h compression/vals16.h

vals16.o : vals16.c compression/vals16.h

bwt.o : bwt.c compression/bwt.h compression/warnmalloc.h

mtf.o : mtf.c compression/mtf.h compression/warnmalloc.h

rle.o : rle.c compression/rle.h

lz77.o : lz77.c compression/lz77.h compression/bwt.h compression/warnmalloc.h

dict.o : dict.c compression/dict.h

huffman.o : huffman.c compression/warnmalloc.h compression/merge_sort.h compression/huffman.h

merge_sort.o : merge_sort.c compression/merge_sort.h compression/warnmalloc.h
