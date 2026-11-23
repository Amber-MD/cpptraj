#!/bin/bash

WORKDIR=`pwd`
BASEDIR=`basename $WORKDIR`

if [ "$BASEDIR" != 'tng' ] ; then
  echo "Execute from inside the 'tng' subdirectory."
  exit 1
fi

SRCDIR=$1
if [ -z "$SRCDIR" ] ; then
  echo "Specify Gromacs TNG source dir."
  exit 1
fi
if [ ! -d "$SRCDIR" ] ; then
  echo "Cannot find TNG source directory $SRCDIR"
  exit 1
fi

UpdateIncludes() {
  if [ -z "$1" ] ; then
    echo "Error: No file given to UpdateIncludes()"
    exit 1
  fi
  sed 's|\.\./\.\./include/||g' $1 > tmp.updateincludes
  mv tmp.updateincludes $1
}

SKIPLIST='CMakeLists.txt README'

for FILE in `git ls-files *` ; do
  echo ""
  BASE=`basename $FILE`
  SKIP=0
  for SKIPFILE in $SKIPLIST ; do
    if [ "$SKIPFILE" = "$BASE" ] ; then
      SKIP=1
      break
    fi
  done
  if [ $SKIP -eq 1 ] ; then
    echo "Skipping $FILE"
    continue
  fi 
  echo "File  : $FILE"
  SRCFILE=''
  for SRC in `find $SRCDIR -name "$BASE"` ; do
    if [ ! -z "$SRCFILE" ] ; then
      echo "Error: Duplicate of $BASE found: $SRC"
      exit 1
    else
      SRCFILE=$SRC
    fi
    if [ -z "$SRCFILE" ] ; then
      echo "Warning: $BASE not found in $SRCDIR"
      continue
    else
      echo "Source: $SRCFILE"
    fi
    M0=`md5sum $FILE | awk '{print $1;}'`
    M1=`md5sum $SRCFILE | awk '{print $1;}'`
    if [ "$M0" != "$M1" ] ; then
      echo "$FILE needs updating."
      echo "cp $SRCFILE => $FILE"
      cp $SRCFILE $FILE
      UpdateIncludes $FILE
      diff -w $FILE $SRCFILE > temp.diff
      cat temp.diff
      rm temp.diff
      #if [ "$BASE" = 'bwlzh.c' ] ; then
      #  echo DEBUG
      #  exit 0
      #fi
    else
      echo "$FILE matches."
    fi
  done
done
