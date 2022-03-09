#!/bin/bash

# Adapt pieces of the Amber cmake build system needed by CPPTRAJ. 

if [ -z "$AMBERHOME" ] ; then
  echo "AMBERHOME is not set."
  exit 1
fi
AMBERCMAKE="$AMBERHOME/cmake"

WORKDIR=`pwd`

# AmberToCpptraj <amber filename> <cpptraj filename>
# This is a 1 to 1 copy of the file, same name, maybe different location.
# If file is already copied, report any differences.
AmberToCpptraj() {
  amberfile=$AMBERCMAKE/$1
  if [ -z "$2" ] ; then
    cpptrajfile=$1
  else
    cpptrajfile=$2
  fi
  echo "AmberToCpptraj: $amberfile $cpptrajfile"
  if [ ! -f "$amberfile" ] ; then
    echo "Amber file $amberfile not present"
    exit 1
  elif [ ! -f "$cpptrajfile" ] ; then
    echo "Copy $amberfile to $cpptrajfile"
  else
    diff $amberfile $cpptrajfile > temp.diff
    if [ -s 'temp.diff' ] ; then
      echo "$amberfile $cpptrajfile diff:"
      cat temp.diff
    fi
    rm temp.diff
  fi
}

# 1 to 1 files
AmberToCpptraj BuildReport.cmake
AmberToCpptraj CheckLinkerFlag.cmake
AmberToCpptraj ColorMessage.cmake
AmberToCpptraj CompilationOptions.cmake
AmberToCpptraj CompilerFlags.cmake
AmberToCpptraj CopyTarget.cmake

exit 0

if ["`basename $WORKDIR`" != 'cmake-cpptraj' ] ; then
  echo "Execute from cmake-cpptraj directory."
  exit 1
fi

# Copy files from Amber that exist here
for FILE in `find . -type f -not -path "*.git/*"` ; do
  NAME=${FILE#./}
  #echo $NAME
  AMBERFILE=$AMBERHOME/cmake/$NAME
  if [ ! -f "$AMBERFILE" ] ; then
    echo "$AMBERFILE not present."
  else
    echo "cp $AMBERFILE -> $NAME"
    #cp $AMBERFILE $NAME
    #git add $NAME
  fi
done
echo ""

# Copy files from Amber that do not exist here
for FILE in `find $AMBERHOME/cmake -type f` ; do
  NAME=${FILE#$AMBERHOME/cmake/}
  #echo $NAME
  if [ ! -f "$NAME" ] ; then
    echo "$NAME is missing."
    #cp $AMBERFILE $NAME
    #git add $NAME
  fi
done
