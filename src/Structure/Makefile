# Structure Makefile
include ../../config.h

include structurefiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(STRUCTURE_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(STRUCTURE_SOURCES) > structuredepend

include structuredepend
