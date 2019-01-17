# Cpptraj Cluster Makefile
include ../../config.h

# Source files to be compiled
SOURCES= \
  Algorithm_HierAgglo.cpp \
  Centroid_Coord.cpp \
  List.cpp \
  Metric_RMS.cpp \
  Node.cpp \
  PairwiseMatrix.cpp \
  PairwiseMatrix_MEM.cpp

# Objects 
OBJECTS=$(SOURCES:.cpp=.o)

# Default target
all: $(OBJECTS)

# Install all targets
install: all

# Dependency targets
../findDepend: ../FindDepend.o 
	$(CXX) -o ../findDepend ../FindDepend.o

depend: ../findDepend
	../findDepend $(SOURCES) > clusterdepend

#dependclean:
#	/bin/rm -f ../FindDepend.o ../findDepend

# Clean/uninstall targets
clean:
	/bin/rm -f $(OBJECTS) 

uninstall: clean

# Header dependencies
include clusterdepend
