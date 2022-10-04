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

# Add source files.
cat > CMakeLists.txt <<EOF
#CMake buildfile for CPPTRAJ $DIR subdirectory.
target_sources(cpptraj_common_obj PRIVATE
EOF

cat > $SOURCEFILE <<EOF
# Files for $DIR subdirectory.
$VARNAME= \\
EOF

LASTFILE=''
for CPPFILE in `ls *.cpp` ; do
  LASTFILE=$CPPFILE
done

for CPPFILE in `ls *.cpp` ; do
  echo "  \${CMAKE_CURRENT_LIST_DIR}/$CPPFILE" >> CMakeLists.txt
  if [ "$CPPFILE" = "$LASTFILE" ] ; then
    echo "  $CPPFILE" >> $SOURCEFILE 
  else
    echo "  $CPPFILE \\" >> $SOURCEFILE
  fi
done

echo ")" >> CMakeLists.txt

# Dependencies
touch $DEPENDFILE

make depend

echo ""
echo "Make sure to add 'include $DIR/$SOURCEFILE' to src/Makefile"
echo "Make sure to add '"$upper"_SOURCEFILES=\$($VARNAME:%.cpp=$DIR/%.cpp)' to src/Makefile"
echo "Make sure to add '\$("$upper"_SOURCEFILES:.cpp=.o)' to the OBJECTS/LIBCPPTRAJ_OBJECTS variables in src/Makefile"
echo "Make sure to add '@echo \$("$upper"_SOURCEFILES)' to the showsources: target in src/Makefile"
echo "Make sure to add '\$("$upper"_SOURCEFILES)' to ./findDepend execution in src/Makefile"
echo "Make sure to add 'cd $DIR && \$(MAKE) clean' to clean: target in src/Makefile"
echo "Make sure to add 'cd $DIR && \$(MAKE) uninstall' to uninstall: target in src/Makefile"
echo "Make sure to add 'add_subdirectory($DIR)' to src/CMakeLists.txt"
echo ""


