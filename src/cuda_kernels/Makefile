# Makefile for cpptraj cuda sources
include ../../external.config.h

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libcpptraj_cuda.a

# General rules

# Source files

CUDA_SOURCES=core_kernels.cu kernel_wrappers.cu GistCudaCalc.cu GistCudaSetup.cu

# Objects

OBJECTS=$(CUDA_SOURCES:.cu=.o)

# Specific rules

all: $(TARGET)

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
	../findDepend $(CUDA_SOURCES) > cudadepend

#dependclean:
#	$(DEL_FILE) FindDepend.o findDepend

include cudadepend
