# CPPTRAJ - Main makefile
# Daniel R. Roe
# 2010-11-18

# Create cpptraj binary in ./src/
all:
	cd src && $(MAKE)

# Create cpptraj binary and move to ./bin/
install_local:
	cd src && $(MAKE) install

# Create cpptraj binary within AmberTools
install:
	cd src && $(MAKE) -f Makefile_at install

# Run Tests
check:
	cd ../../test/cpptraj/ && $(MAKE) test

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

