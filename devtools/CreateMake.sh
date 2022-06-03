#!/bin/bash

# Create files needed for build systems in /src subdirectories.

WORKDIR=`pwd`

DIRNAME=`dirname $WORKDIR`

BASENAME=`basename $DIRNAME`

if [ "$BASENAME" != 'src' ] ; then
  echo "Error: Must be executed in a subdirectory of the CPPTRAJ src directory."
  exit 1
fi


