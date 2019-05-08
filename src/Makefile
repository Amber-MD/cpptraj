# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

DEL_FILE = /bin/rm -f

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o)

# Default target: cpptraj only
all: cpptraj$(SFX)$(EXE)

# Install all targets
install: $(INSTALL_TARGETS) 

# Cpptraj ------------------------------
install_cpptraj: cpptraj$(SFX)$(EXE)
	mv cpptraj$(SFX)$(EXE) $(CPPTRAJBIN)

cpptraj$(SFX)$(EXE): $(OBJECTS) $(FFT_TARGET) $(READLINE_TARGET) $(CUDA_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET)
	$(CXX) -o cpptraj$(SFX)$(EXE) $(OBJECTS) $(CUDA_TARGET) $(FFT_TARGET) $(READLINE_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

# libcpptraj ---------------------------
# Rule to make libcpptraj-specific objects
%.LIBCPPTRAJ.o : %.cpp
	$(CXX) $(DIRECTIVES) -DLIBCPPTRAJ $(INCLUDE) $(CXXFLAGS) -c -o $@ $<

libcpptraj: $(LIBCPPTRAJ_TARGET)

$(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX): $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET)
	$(CXX) -shared -o $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX) $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET) $(CPPTRAJ_LIB) $(LDFLAGS)

$(CPPTRAJLIB)/libcpptraj.a: $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET)
	ar rv $(CPPTRAJLIB)/libcpptraj.a $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET)

libstatic: $(CPPTRAJLIB)/libcpptraj.a

nolibcpptraj:
	@echo ""
	@echo "Error: Cannot build libcpptraj; re-configure with '-shared'"
	@echo ""
	@exit 1

# nproc --------------------------------
nproc:
	cd ../util && $(MAKE)

# Internal object/library targets
pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

readline/libreadline.a: ../external.config.h
	cd readline && $(MAKE) all

xdrfile/libxdrfile.a: ../external.config.h
	cd xdrfile && $(MAKE) all

arpack/libarpack.a: ../external.config.h
	cd arpack && $(MAKE) install

noreadline:
	@echo ""
	@echo "Skipping bundled READLINE build"
	@echo ""

noxdrfile:
	@echo ""
	@echo "Skipping bundled XDRFILE build"
	@echo ""

noarpack:
	@echo ""
	@echo "Skipping bundled ARPACK build"
	@echo ""

cuda_kernels/libcpptraj_cuda.a: ../external.config.h
	cd cuda_kernels && $(MAKE) all

# Dependency targets
findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) > cpptrajdepend

dependclean:
	$(DEL_FILE) FindDepend.o findDepend

# Clean/uninstall targets
clean:
	$(DEL_FILE) $(OBJECTS) $(FFT_TARGET) cpptraj$(SFX) *.LIBCPPTRAJ.o libcpptraj$(SHARED_SUFFIX)
	cd readline && $(MAKE) clean
	cd xdrfile && $(MAKE) clean
	cd arpack && $(MAKE) clean
	cd cuda_kernels && $(MAKE) clean

uninstall:
	$(DEL_FILE) $(CPPTRAJBIN)/cpptraj$(SFX)$(EXE)
	$(DEL_FILE) $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX)
	cd readline && make uninstall
	cd xdrfile && make uninstall
	cd cuda_kernels && make uninstall

# Header dependencies
include cpptrajdepend
