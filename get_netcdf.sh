#!/bin/bash

WORKDIR=`pwd`
# Attempt to download and install a copy of NetCDF
if [ -z "$SRCTAR" -o -z "$SRCDIR" -o -z "$URL" ] ; then
  echo "Error: Script variables are empty."
  exit 1
fi

# Get netcdf if necessary
if [ ! -f "$SRCTAR" ] ; then
  WGET=`which wget`

  if [ -z "$WGET" ] ; then
    echo "Error: 'wget' not found. Cannot download NetCDF"
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
  tar -zxf $SRCTAR
  if [ $? -ne 0 -o ! -d "$SRCDIR" ] ; then
    echo "Error: Could not unpack $SRCTAR"
    exit 1
  fi
fi
cd $SRCDIR

# Configure
echo ""
echo "CC=$CC"
echo "CFLAGS=$CFLAGS"
echo "PREFIX=$PREFIX"
echo -n "Configuring NetCDF... "
./configure CC="$CC" CFLAGS="$CFLAGS" \
  --prefix=$PREFIX --disable-netcdf-4 --disable-dap $windows_hostflag \
  --disable-shared --disable-doxygen > ../netcdf_config.log 2>&1
if [ $? -ne 0 ] ; then
  echo "Failed."
  echo "Check $WORKDIR/netcdf_config.log for errors."
  exit 1
else
  echo "Success."
fi

# Build
echo -n "Compiling NetCDF... "
make > ../netcdf_compile.log 2>&1
if [ $? -ne 0 ] ; then
  echo "Build failed."
  echo "Check $WORKDIR/netcdf_compile.log for errors."
  exit 1
fi

# Install
make install >> ../netcdf_compile.log 2>&1
if [ $? -ne 0 ] ; then
  echo "Install failed."
  echo "Check $WORKDIR/netcdf_compile.log for errors."
  exit 1
fi
echo "Success."

exit 0
