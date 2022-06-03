#!/bin/bash

# Create files needed for build systems in /src subdirectories.

WORKDIR=`pwd`

DIRNAME=`dirname $WORKDIR`

BASENAME=`basename $DIRNAME`

if [ "$BASENAME" != 'src' ] ; then
  echo "Error: Must be executed in a subdirectory of the CPPTRAJ src directory."
  exit 1
fi

DIR=`basename $WORKDIR`
name=`echo "$DIR" | awk '{print tolower($0)}'`
upper=`echo "$DIR" | awk '{print toupper($0)}'`
echo "Dir= $DIR, name= $name"

SOURCEFILE="$name"files
DEPENDFILE="$name"depend
VARNAME="$upper"_SOURCES

echo "File containing Makefile sources: $SOURCEFILE"
echo "File containing Makefile depends: $DEPENDFILE"

#if [ ! -f 'Makefile' ] ; then
  # Create the Makefile
  cat > Makefile <<EOF
# $DIR Makefile
include ../../config.h

include $SOURCEFILE

DEL_FILE = /bin/rm -f

# Objects
OBJECTS=\$($VARNAME:.cpp=.o)

# Default target: objects
all: \$(OBJECTS)

clean:
	\$(DEL_FILE) *.o

uninstall: clean

# Dependency targets
../findDepend:
	cd ../ && \$(MAKE) findDepend

depend: ../findDepend
	../findDepend \$($VARNAME) > $DEPENDFILE

include $DEPENDFILE
EOF

echo "Make sure to add 'include $DIR/$SOURCEFILE' to src/Makefile"
