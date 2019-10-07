# Makefile for TNG bundled with cpptraj
include ../../external.config.h

CFLAGS += -I .

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libtngfile.a

# Source files
SOURCES=md5.c tng_io.c tng_compress.c coder.c warnmalloc.c fixpoint.c

# Objects
OBJECTS=$(SOURCES:.c=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

test: $(TARGET) tng_io_testing.o
	$(CC) -o a.out tng_io_testing.c $(TARGET) -lm $(ZLIB_FLAG)
	./a.out

clean:
	$(DEL_FILE) *.o $(TARGET) tng_io_testing.o a.out

uninstall: clean

# Dependencies
md5.o : md5.c tng/md5.h

tng_io.o : tng_io.c tng/tng_io.h tng/tng_io_fwd.h compression/tng_compress.h tng/version.h

tng_compress.o : tng_compress.c compression/tng_compress.h compression/coder.h compression/fixpoint.h

coder.o : coder.c compression/coder.h compression/tng_compress.h compression/bwlzh.h compression/warnmalloc.h

warnmalloc.o : warnmalloc.c compression/warnmalloc.h compression/tng_compress.h

fixpoint.o : fixpoint.c compression/fixpoint.h 
