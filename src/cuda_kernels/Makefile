# Makefile for cpptraj cuda sources
include ../../external.config.h

# Variables
DEL_FILE      = /bin/rm -f
AR            = ar cqs
TARGET        = libcpptraj_cuda.a

# General rules

# Source files

CUDA_SOURCES=core_kernels.cu kernel_wrappers.cu

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

# Dependencies
core_kernels.o: core_kernels.cu

kernel_wrappers.o: kernel_wrappers.cu
