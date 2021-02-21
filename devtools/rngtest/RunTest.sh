#!/bin/bash

COUNT=100000000

make rngtest

for CPPTRAJTYPE in 0 1 2 ; do
  echo "Cpptraj type $CPPTRAJTYPE"
  TMPFILE=temp.$CPPTRAJTYPE
  ./rngtest -t $COUNT -S 1 -r $CPPTRAJTYPE > $TMPFILE
  for TEST in 0 1 2 3 4 ; do
    dieharder -d $TEST -f $TMPFILE -g 202
  done
  echo ""
done

