# Energy Makefile
include ../../config.h

include energyfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(ENERGY_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(ENERGY_SOURCES) > energydepend

include energydepend
