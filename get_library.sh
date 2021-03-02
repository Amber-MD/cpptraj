#!/bin/bash

WORKDIR=`pwd`

# Attempt to download and install a copy of library 
if [ -z "$SRCTAR" -o -z "$URL" -o -z "$LIBNAME" ] ; then
  echo "Error: Script variables are empty."
  exit 1
fi

# First ask if we want to get the library
echo "Should CPPTRAJ attempt to build its own $LIBNAME? {y|n}: "
read yesno
while [ "$yesno" != 'y' -a "$yesno" != 'n' ] ; do
  echo "    Please enter 'y' or 'n': "
  read yesno
done
if [ "$yesno" = 'n' ] ; then
  exit 1
fi

# Get configure options
CONFIGOPTS=''
while [ ! -z "$1" ] ; do
  CONFIGOPTS="$CONFIGOPTS $1"
  shift
done

echo "    CONFIGOPTS: $CONFIGOPTS"

CONFIGURE_LOG=$WORKDIR/$LIBNAME"_config.log"
COMPILE_LOG=$WORKDIR/$LIBNAME"_compile.log"

# Get library if necessary
if [ ! -f "$SRCTAR" ] ; then
  echo "    Downloading $LIBNAME..."
  WGET=`which wget`

  if [ -z "$WGET" ] ; then
    echo "Error: 'wget' not found. Cannot download $LIBNAME"
    exit 1
  fi

  $WGET $URL
  if [ $? -ne 0 -o ! -f "$SRCTAR" ] ; then
    echo "Error: Could not download $URL"
    exit 1
  fi
fi

# Get SRCDIR if necessary
if [ -z "$SRCDIR" ] ; then
  FIRSTFILE=`tar -tzf $SRCTAR | head -1`
  if [ ! -z "$FIRSTFILE" ] ; then
    SRCDIR=`dirname $FIRSTFILE`
  fi
fi

# Unpack
if [ -z "$SRCDIR" ] ; then
  echo "Error: SRCDIR is empty."
  exit 1
fi
if [ ! -d "$SRCDIR" ] ; then
  echo "    Unpacking $LIBNAME..."
  tar -zxf $SRCTAR
  if [ $? -ne 0 -o ! -d "$SRCDIR" ] ; then
    echo "Error: Could not unpack $SRCTAR"
    exit 1
  fi
fi
cd $SRCDIR

# Configure
#echo ""
#echo "CC=$CC"
#echo "CFLAGS=$CFLAGS"
#echo "PREFIX=$PREFIX"
echo -n "    Configuring $LIBNAME... "
if [ -f 'configure' ] ; then
  # Run configure
  ./configure CC="$CC" CFLAGS="$CFLAGS" --prefix=$PREFIX $CONFIGOPTS > $CONFIGURE_LOG 2>&1
  if [ $? -ne 0 ] ; then
    echo "Failed."
    echo "Check $CONFIGURE_LOG for errors."
    exit 1
  fi
elif [ -f 'Makefile' ] ; then
  # No configure - try to modify Makefile on the fly.
  if [ ! -f 'Makefile.original' ] ; then
    cp Makefile Makefile.original
  fi
  awk -v cc="$CC" -v cflags="$CFLAGS" -v prefix="$PREFIX" '{
    if (index($1,"CC=")!=0)
      printf("CC=%s\n", cc);
    else if (index($1,"CFLAGS=")!=0)
      printf("CFLAGS=%s\n", cflags);
    else if (index($1,"PREFIX=")!=0)
      printf("PREFIX=%s\n", prefix);
    else
      print $0;
  }' Makefile.original > Makefile
  if [ $? -ne 0 ] ; then
    exit 1
  fi
else
  echo "Error: $LIBNAME has no 'configure' or 'Makefile'."
  exit 1
fi
echo "Success."

# Build
echo -n "    Compiling $LIBNAME... "
make clean > $COMPILE_LOG 2>&1
make > $COMPILE_LOG 2>&1
if [ $? -ne 0 ] ; then
  echo "Build failed."
  echo "Check $COMPILE_LOG for errors."
  exit 1
fi

# Install
make install >> $COMPILE_LOG 2>&1
if [ $? -ne 0 ] ; then
  echo "Install failed."
  echo "Check $COMPILE_LOG for errors."
  exit 1
fi
echo "Success."

exit 0
