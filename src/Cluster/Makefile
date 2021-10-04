# Cluster code library Makefile
include ../../config.h

include clusterfiles

DEL_FILE = /bin/rm -f
AR            = ar cqs
TARGET        = libcpptraj_cluster.a

# Objects
OBJECTS=$(CLUSTER_SOURCES:.cpp=.o)

#CLUSTER_TARGET=Cluster/libcpptraj_cluster.a

# Default target: library
all: $(TARGET) 

# Build library
$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

clean:
	$(DEL_FILE) *.o $(TARGET)

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && $(MAKE) findDepend

depend: ../findDepend
	../findDepend $(CLUSTER_SOURCES) > clusterdepend

#dependclean:
#	$(DEL_FILE) FindDepend.o findDepend

include clusterdepend
