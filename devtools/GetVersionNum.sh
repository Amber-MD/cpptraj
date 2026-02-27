#!/bin/bash

EXEDIR=`dirname $0`

if [ ! -d "$EXEDIR" ] ; then
  echo "Error: Could not determine execution directory."
  exit 1
fi

cd $EXEDIR

VFILE=../src/Version.h
if [ ! -f "$VFILE" ] ; then
  echo "Version file $VFILE not found."
  exit 1
fi

VERSION=`awk '{if ($2 == "CPPTRAJ_INTERNAL_VERSION") {gsub(/\"/,""); printf ("%s", $3);}}' $VFILE`
echo $VERSION
exit 0
