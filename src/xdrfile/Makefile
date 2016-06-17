include ../../config.h

# Variables

AR            = ar cqs
DEL_FILE      = /bin/rm -f
TARGET        = libxdrfile.a

# General rules

#%.o : %.cu
#	$(NVCC) $(DBGFLAGS) -c $(NVCCFLAGS) -o $@ $<

# Source files

XDR_SOURCES=xdrfile.c xdrfile_trr.c xdrfile_xtc.c

# Objects

OBJECTS=$(XDR_SOURCES:.c=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

clean:
	/bin/rm -f *.o $(TARGET)

# Dependencies
xdrfile.o: xdrfile.c xdrfile.h

xdrfile_trr.o: xdrfile_trr.c xdrfile_trr.h xdrfile.h

xdrfile_xtc.o: xdrfile_xtc.c xdrfile_xtc.h
