# CPPTRAJ - Main Makefile
# Daniel R. Roe
# 2010-11-18
# Revised 2015-03-09

# Create standalone cpptraj binary in ./src/
all:
	cd src && $(MAKE)

# Create standalone cpptraj binary
install_local:
	cd src && $(MAKE) install

# Create cpptraj/ambpdb binaries within AmberTools
serial:
	cd src && $(MAKE) -f Makefile_at install

# Create OpenMP cpptraj binary within AmberTools
openmp:
	cd src && $(MAKE) -f Makefile_at install_openmp

# Create MPI cpptraj binary within AmberTools
parallel:
	cd src && $(MAKE) -f Makefile_at install_mpi

# Create libcpptraj within AmberTools 
libcpptraj:
	cd src && $(MAKE) -f Makefile_at libcpptraj

# Run Tests
check:
	cd ../../test/cpptraj/ && $(MAKE) test

check_local:
	cd ../../test/cpptraj/ && $(MAKE) test.standalone

docs: src/cpptraj.Doxyfile
	cd src && doxygen cpptraj.Doxyfile

# Clean up for standalone
clean_local:
	cd src && $(MAKE) clean

# Clean up for AmberTools
clean:
	cd src && $(MAKE) -f Makefile_at clean

# Uninstall for AmberTools 
uninstall:
	cd src && $(MAKE) -f Makefile_at uninstall
