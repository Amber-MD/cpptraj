# Cluster Makefile
include ../../config.h

include clusterfiles

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=$(CLUSTER_SOURCES:.cpp=.o)

# Default target: objects
all: $(OBJECTS)

clean:
	$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(CLUSTER_SOURCES) > clusterdepend

include clusterdepend
