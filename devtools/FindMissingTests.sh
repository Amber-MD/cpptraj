#!/bin/bash

if [ ! -f 'Makefile' ] ; then
  echo "Should be executed in directory where main test Makefile is."
  exit 1
fi

for DIR in `ls -d Test_*` ; do
  if [ -d "$DIR" ] ; then
    IN_MAKEFILE=`grep $DIR Makefile`
    if [ -z "$IN_MAKEFILE" ] ; then
      echo "$DIR not in Makefile."
    fi
  fi
done
