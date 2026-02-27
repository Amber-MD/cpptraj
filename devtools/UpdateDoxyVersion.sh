#!/bin/bash

# Update the cpptraj.Doxyfile version number

EXEDIR=`dirname $0`

if [ ! -d "$EXEDIR" ] ; then
  echo "Error: Could not determine execution directory."
  exit 1
fi

cd $EXEDIR
DFILE=../src/cpptraj.Doxyfile
if [ ! -f "$DFILE" ] ; then
  echo "Doxygen file $DFILE not found."
  exit 1
fi

VERSION=`./GetVersionNum.sh`
if [ $? -ne 0 ] ; then
  echo "Error: Could not get version number."
  exit 1
fi

if [ -f 'tmp.cpptraj.Doxyfile' ] ; then
  rm tmp.cpptraj.Doxyfile
fi
awk -v vstring="$VERSION" '{
  if ($1 == "PROJECT_NUMBER")
    printf("PROJECT_NUMBER         = %s\n", vstring);
  else
    print $0;
}' $DFILE > tmp.cpptraj.Doxyfile
mv tmp.cpptraj.Doxyfile $DFILE

exit 0
