# Makefile for xdrfile bundled with cpptraj
include xdrfile.config.h

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libxdrfile.a

# Source files
XDR_SOURCES=xdrfile.c xdr_seek.c xdrfile_trr.c xdrfile_xtc.c

# Objects
OBJECTS=$(XDR_SOURCES:.c=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

test: $(TARGET) xdrfile_c_test.o
	$(CC) -o a.out xdrfile_c_test.c $(TARGET) -lm
	./a.out

clean:
	$(DEL_FILE) *.o $(TARGET) a.out test.trr test.xtc test.xdr

uninstall: clean

# Dependencies
xdrfile.o: xdrfile.c xdrfile.h

xdrfile_trr.o: xdrfile_trr.c xdrfile_trr.h xdrfile.h

xdrfile_xtc.o: xdrfile_xtc.c xdrfile_xtc.h xdrfile.h

xdr_seek.o: xdr_seek.c xdr_seek.h xdrfile.h

xdrfile_c_test.o: xdrfile_c_test.c xdrfile_xtc.h xdrfile_trr.h xdrfile.h
