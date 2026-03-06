# Parm Makefile
include ../../config.h

include parmfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(PARM_SOURCES:.cpp=.o)

# General rules
.cpp.o:
	$(VB)echo CXX $<
	$(VB)$(CXX) $(DIRECTIVES) $(CPPTRAJ_INC) $(CXXFLAGS) -c -o $@ $<

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(PARM_SOURCES) > parmdepend

include parmdepend
