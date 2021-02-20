#!/bin/bash

COUNT=1000000

make rngtest

for CPPTRAJTYPE in 0 1 ; do
  echo "Cpptraj type $CPPTRAJTYPE"
  TMPFILE=temp.$CPPTRAJTYPE
  ./rngtest -t $COUNT -S 1 -r $CPPTRAJTYPE > $TMPFILE
  dieharder -d 0 -f $TMPFILE -g 202
  echo ""
done

