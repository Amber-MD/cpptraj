include ../../config.h

# Variables

AR            = ar cqs
DEL_FILE      = /bin/rm -f
TARGET        = libcpptraj_cuda.a

# General rules

%.o : %.cu
	$(NVCC) $(DBGFLAGS) -c $(NVCCFLAGS) -o $@ $<

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
	/bin/rm -f *.o $(TARGET)

# Dependencies
core_kernels.o: core_kernels.cu

kernel_wrappers.o: kernel_wrappers.cu
