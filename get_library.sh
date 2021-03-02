#!/bin/bash

WORKDIR=`pwd`
# Attempt to download and install a copy of library 
if [ -z "$SRCTAR" -o -z "$SRCDIR" -o -z "$URL" -o -z "$LIBNAME" ] ; then
  echo "Error: Script variables are empty."
  exit 1
fi

# Get configure options
CONFIGOPTS=''
while [ ! -z "$1" ] ; do
  CONFIGOPTS="$CONFIGOPTS $1"
  shift
done

echo "CONFIGOPTS: $CONFIGOPTS"

CONFIGURE_LOG=$WORKDIR/$LIBNAME"_config.log"
COMPILE_LOG=$WORKDIR/$LIBNAME"_compile.log"

# Get library if necessary
if [ ! -f "$SRCTAR" ] ; then
  echo "Downloading $LIBNAME..."
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

# Unpack
if [ ! -d "$SRCDIR" ] ; then
  echo "Unpacking $LIBNAME..."
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
echo -n "Configuring $LIBNAME... "
./configure CC="$CC" CFLAGS="$CFLAGS" --prefix=$PREFIX $CONFIGOPTS > $CONFIGURE_LOG 2>&1
if [ $? -ne 0 ] ; then
  echo "Failed."
  echo "Check $CONFIGURE_LOG for errors."
  exit 1
else
  echo "Success."
fi

# Build
echo -n "Compiling $LIBNAME... "
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
