# CPPTRAJ - Main makefile
# Daniel R. Roe
# 2010-11-18

# Create cpptraj binary in ./src/
all:
	cd src && $(MAKE)

# Create cpptraj binary and move to ./bin/
install_local:
	cd src && $(MAKE) install

# Create cpptraj/ambpdb binaries within AmberTools
install:
	cd src && $(MAKE) -f Makefile_at install

# Create OpenMP cpptraj binary within AmberTools
install_openmp:
	cd src && $(MAKE) -f Makefile_at install_openmp

# Create MPI cpptraj binary within AmberTools
install_mpi:
	cd src && $(MAKE) -f Makefile_at install_mpi

# Run Tests
check:
	cd ../../test/cpptraj/ && $(MAKE) test

check_local:
	cd ../../test/cpptraj/ && $(MAKE) test.standalone

# Clean up
clean_local:
	cd src && $(MAKE) clean

docs: src/cpptraj.Doxyfile
	cd src && doxygen cpptraj.Doxyfile

# Clean up for AmberTools
clean:
	cd src && $(MAKE) -f Makefile_at clean

# Clean up lapack, arpack, and blas libraries for standalone
libclean:
	cd src && $(MAKE) libclean

# Remove cpptraj binary from $AMBERHOME/bin
uninstall:
	cd src && $(MAKE) -f Makefile_at uninstall

