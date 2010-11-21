# CPPTRAJ - Main makefile
# Daniel R. Roe
# 2010-11-18

# Create cpptraj binary in ./src/
all:
	cd src && $(MAKE)

# Create cpptraj binary and move to ./bin/
install:
	cd src && $(MAKE) install

# Run Tests
check:
	cd test && $(MAKE) test

# Clean up
clean:
	cd src && $(MAKE) clean
