#!/bin/bash

# Attempt to download and install a copy of NetCDF
URL='ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.6.1.tar.gz'

# Necessary programs
WGET=`which wget`

if [ -z "$WGET" ] ; then
  echo "Error: 'wget' not found. Cannot download NetCDF"
  exit 1
fi


