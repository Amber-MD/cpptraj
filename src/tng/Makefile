# Makefile for TNG bundled with cpptraj
include ../../external.config.h

CFLAGS += -I .

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libtngfile.a

# Source files
SOURCES=md5.c tng_io.c

# Objects
OBJECTS=$(SOURCES:.c=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

test: $(TARGET) tng_io_testing.o
	$(CC) -o a.out tng_io_testing.c $(TARGET) -lm
	./a.out

clean:
	$(DEL_FILE) *.o $(TARGET) tng_io_testing.o a.out

uninstall: clean

# Dependencies
md5.o : md5.c tng/md5.h

tng_io.o : tng_io.c tng/tng_io.h tng/tng_io_fwd.h compression/tng_compress.h tng/version.h
