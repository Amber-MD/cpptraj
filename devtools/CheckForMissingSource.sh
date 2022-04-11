#!/bin/bash

MAKE_SOURCES=`ls *files`
if [ -z "$MAKE_SOURCES" ] ; then
  MAKE_SOURCES=Makefile
fi
CMAKE_SOURCES=CMakeLists.txt

echo "Make sources  : $MAKE_SOURCES"
echo "Cmake sources : $CMAKE_SOURCES"

if [ ! -f "$MAKE_SOURCES" ] ; then
  echo "Make sources not found."
  exit 1
fi

if [ ! -f "$CMAKE_SOURCES" ] ; then
  echo "Cmake sources not found."
  exit 1
fi

SOURCES=`ls *.cpp *.c *.F90 2> /dev/null`

for FILE1 in $MAKE_SOURCES $CMAKE_SOURCES ; do
  for FILE2 in $SOURCES ; do
    if [ -z "`grep $FILE2 $FILE1`" ] ; then
      echo "$FILE2 appears to be missing from $FILE1"
    fi
  done
done
