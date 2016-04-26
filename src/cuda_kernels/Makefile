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

SOURCES=setup_routines.cpp

# Objects

OBJECTS=$(CUDA_SOURCES:.cu=.o) $(SOURCES:.cpp=.o)

# Specific rules

all: $(TARGET)

$(TARGET): $(OBJECTS)
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)

clean:
	/bin/rm -f *.o $(TARGET)

# Dependencies
setup_routines.o: setup_routines.cpp ../Action_Closest.h

core_kernels.o: core_kernels.cu

kernel_wrappers.o: kernel_wrappers.cu
