# Cpptraj Cluster Makefile
include ../../config.h

# Variables
AR       = ar cqs
DEL_FILE = /bin/rm -f
TARGET   = libcpptraj_cluster.a

# Source files to be compiled
SOURCES= \
  Algorithm_HierAgglo.cpp \
  Centroid_Coord.cpp \
  Control.cpp \
  List.cpp \
  Metric_RMS.cpp \
  Node.cpp \
  PairwiseMatrix.cpp \
  PairwiseMatrix_MEM.cpp

# Objects 
OBJECTS=$(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Library target
$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

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
	$(DEL_FILE) $(OBJECTS)

uninstall: clean

# Header dependencies
include clusterdepend
