# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o)

include ambpdbfiles

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o)

CUDA_SOURCES=cuda_kernels/core_kernels.cu cuda_kernels/kernel_wrappers.cu

# Default target: cpptraj only
all: cpptraj$(SFX)$(EXE)

# Install all targets
install: $(INSTALL_TARGETS) 

# Cpptraj ------------------------------
install_cpptraj: cpptraj$(SFX)$(EXE)
	mv cpptraj$(SFX)$(EXE) $(CPPTRAJBIN)

cpptraj$(SFX)$(EXE): $(OBJECTS) $(FFT_DEPEND) $(READLINE_TARGET) $(CUDA_TARGET) $(XDRFILE_TARGET)
	$(CXX) -o cpptraj$(SFX)$(EXE) $(OBJECTS) $(READLINE) $(CUDA_TARGET) $(FFT_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

# AmbPDB -------------------------------
install_ambpdb: ambpdb$(EXE)
	mv ambpdb$(EXE) $(CPPTRAJBIN)

ambpdb$(EXE): $(AMBPDB_OBJECTS) $(XDRFILE_TARGET)
	$(CXX) -o ambpdb$(EXE) $(AMBPDB_OBJECTS) $(LDFLAGS)

# libcpptraj ---------------------------
# Rule to make libcpptraj-specific objects
%.LIBCPPTRAJ.o : %.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) -DLIBCPPTRAJ $(DBGFLAGS) $<

libcpptraj: $(LIBCPPTRAJ_TARGET)

$(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX): $(LIBCPPTRAJ_OBJECTS) $(FFT_DEPEND)
	$(CXX) -shared -o $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX) $(LIBCPPTRAJ_OBJECTS) $(FFT_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

nolibcpptraj:
	@echo ""
	@echo "Error: Cannot build libcpptraj; re-configure with '-shared'"
	@echo ""
	@exit 1

# nproc --------------------------------
nproc:
	cd ../util && $(MAKE)

# Additional library targets
readline/libreadline.a: $(READLINE_HOME)/Makefile
	cd $(READLINE_HOME) && $(MAKE) all

xdrfile/libxdrfile.a: $(XDRFILE_HOME)/Makefile
	cd $(XDRFILE_HOME) && $(MAKE) all

noreadline:
	@echo ""
	@echo "Skipping bundled READLINE build"
	@echo ""

noxdrfile:
	@echo ""
	@echo "Skipping bundled XDRFILE build"
	@echo ""

cuda_kernels/libcpptraj_cuda.a: $(CUDA_SOURCES)
	cd cuda_kernels && $(MAKE) all

# Individual file targets
$(FFT_DEPEND):  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

ReadLine.o: ReadLine.cpp
	$(CXX) $(CXXFLAGS) -I$(READLINE_HOME) -c -o $@ ReadLine.cpp

# Dependency targets
findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) AmbPDB.cpp > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

# Clean/uninstall targets
clean:
	/bin/rm -f $(OBJECTS) $(FFT_DEPEND) cpptraj$(SFX) AmbPDB.o ambpdb *.LIBCPPTRAJ.o libcpptraj$(SHARED_SUFFIX)
	cd readline && $(MAKE) clean
	cd xdrfile && $(MAKE) clean
	cd cuda_kernels && $(MAKE) clean

uninstall:
	/bin/rm -f $(CPPTRAJBIN)/cpptraj $(CPPTRAJBIN)/cpp_ambpdb readline/libreadline.a $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX) cuda_kernels/libcpptraj_cuda.a

# Header dependencies
include cpptrajdepend
