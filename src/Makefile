# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

DEL_FILE = /bin/rm -f

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o)

CLUSTER_TARGET=Cluster/libcpptraj_cluster.a

# Default target: cpptraj only
all: cpptraj$(SFX)$(EXE)

# Install all targets
install: $(INSTALL_TARGETS)

# Cpptraj ------------------------------
install_cpptraj: cpptraj$(SFX)$(EXE)
	mv cpptraj$(SFX)$(EXE) $(CPPTRAJBIN)

cpptraj$(SFX)$(EXE): $(OBJECTS) $(FFT_TARGET) $(READLINE_TARGET) $(CUDA_TARGET) $(CLUSTER_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET) $(TNGFILE_TARGET)
	$(CXX) -o cpptraj$(SFX)$(EXE) $(OBJECTS) $(CUDA_TARGET) $(CLUSTER_TARGET) $(FFT_TARGET) $(READLINE_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

# libcpptraj ---------------------------
# Rule to make libcpptraj-specific objects
%.LIBCPPTRAJ.o : %.cpp
	$(CXX) $(DIRECTIVES) -DLIBCPPTRAJ $(INCLUDE) $(CXXFLAGS) -c -o $@ $<

libcpptraj: $(LIBCPPTRAJ_TARGET)

$(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX): $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET) $(CLUSTER_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET) $(TNGFILE_TARGET)
	$(CXX) -shared -o $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX) $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET) $(CLUSTER_TARGET) $(CPPTRAJ_LIB) $(LDFLAGS)

# Data directory -----------------------
install_dat:
	cp ../dat/ATOMTYPE_CHECK.TAB $(CPPTRAJDAT)/
	cp ../dat/BONDTYPE_CHECK.TAB $(CPPTRAJDAT)/
	cp ../dat/Karplus.txt        $(CPPTRAJDAT)/
	cp ../dat/README             $(CPPTRAJDAT)/

# Static libraries ---------------------
#$(CPPTRAJLIB)/libcpptraj.a: $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET)
#	ar rv $(CPPTRAJLIB)/libcpptraj.a $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CUDA_TARGET)

$(CPPTRAJLIB)/libcpptraj_core.a: $(LIBCPPTRAJ_CORE_OBJECTS)
	ar rv $(CPPTRAJLIB)/libcpptraj_core.a $(LIBCPPTRAJ_CORE_OBJECTS)

$(CPPTRAJLIB)/libcpptraj_file.a: $(LIBCPPTRAJ_FILE_OBJECTS)
	ar rv $(CPPTRAJLIB)/libcpptraj_file.a $(LIBCPPTRAJ_FILE_OBJECTS)

$(CPPTRAJLIB)/libcpptraj_traj.a: $(LIBCPPTRAJ_TRAJ_OBJECTS) $(XDRFILE_OBJECTS) $(TNGFILE_OBJECTS) $(VMDPLUGIN_OBJECTS)
	ar rv $(CPPTRAJLIB)/libcpptraj_traj.a $(LIBCPPTRAJ_TRAJ_OBJECTS) $(XDRFILE_OBJECTS) $(TNGFILE_OBJECTS) $(VMDPLUGIN_OBJECTS)

$(CPPTRAJLIB)/libcpptraj_parm.a: $(LIBCPPTRAJ_PARM_OBJECTS)
	ar rv $(CPPTRAJLIB)/libcpptraj_parm.a $(LIBCPPTRAJ_PARM_OBJECTS)

libstatic: $(CPPTRAJLIB)/libcpptraj_core.a $(CPPTRAJLIB)/libcpptraj_traj.a \
           $(CPPTRAJLIB)/libcpptraj_file.a $(CPPTRAJLIB)/libcpptraj_parm.a

nolibcpptraj:
	@echo ""
	@echo "Error: Cannot build libcpptraj; re-configure with '-shared'"
	@echo ""
	@exit 1

# Internal object/library targets
pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

readline/libreadline.a: ../external.config.h
	cd readline && $(MAKE) all

xdrfile/libxdrfile.a: ../external.config.h
	cd xdrfile && $(MAKE) all

arpack/libarpack.a: ../external.config.h
	cd arpack && $(MAKE) install

tng/libtng_io.a: ../external.config.h
	cd tng && $(MAKE) all

noreadline:
	@echo ""
	@echo "Skipping bundled READLINE build"
	@echo ""

noxdrfile:
	@echo ""
	@echo "Skipping bundled XDRFILE build"
	@echo ""

notngfile:
	@echo ""
	@echo "Skipping bundled TNG build"
	@echo ""

noarpack:
	@echo ""
	@echo "Skipping bundled ARPACK build"
	@echo ""

cuda_kernels/libcpptraj_cuda.a: ../external.config.h
	cd cuda_kernels && $(MAKE) all

$(CLUSTER_TARGET)::
	cd Cluster && $(MAKE) all

# Dependency targets
findDepend: FindDepend.cpp FindDepend.o 
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
	cd Cluster && $(MAKE) clean
	cd tng && $(MAKE) clean

uninstall_lib:
	$(DEL_FILE) $(CPPTRAJLIB)/*.a
	$(DEL_FILE) $(CPPTRAJLIB)/*.la
	$(DEL_FILE) $(CPPTRAJLIB)/*.so*
	$(DEL_FILE) $(CPPTRAJLIB)/*.settings
	$(DEL_FILE) -r $(CPPTRAJLIB)/cmake
	$(DEL_FILE) -r $(CPPTRAJLIB)/pkgconfig

uninstall_inc:
	$(DEL_FILE) $(CPPTRAJINC)/*.h
	$(DEL_FILE) $(CPPTRAJINC)/*.f
	$(DEL_FILE) $(CPPTRAJINC)/*.f03

uninstall: uninstall_lib uninstall_inc
	$(DEL_FILE) $(CPPTRAJBIN)/cpptraj$(SFX)$(EXE)
	$(DEL_FILE) $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX)
	cd readline && $(MAKE) uninstall
	cd xdrfile && $(MAKE) uninstall
	cd cuda_kernels && $(MAKE) uninstall
	cd Cluster && $(MAKE) uninstall

# Header dependencies
include cpptrajdepend
