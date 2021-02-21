#!/bin/bash

COUNT=100000000

make rngtest

for CPPTRAJTYPE in 0 1 2 3 ; do
  OUTFILE=results.$CPPTRAJTYPE
  echo "Cpptraj type $CPPTRAJTYPE"
  echo "Cpptraj type $CPPTRAJTYPE" > $OUTFILE
  TMPFILE=temp.$CPPTRAJTYPE
  ./rngtest -t $COUNT -S 1 -r $CPPTRAJTYPE > $TMPFILE
  for TEST in 0 1 2 3 4 ; do
    echo "Test $TEST"
    dieharder -d $TEST -f $TMPFILE -g 202 >> $OUTFILE
  done
  echo ""
done

echo "Diehard numbers"
echo "Diehard numbers" > results.diehard
dieharder -S 1 -B -o -t $COUNT > diehard.numbers
if [ -f 'diehard.numbers' ] ; then
  for TEST in 0 1 2 3 4 ; do
    echo "Test $TEST"
    dieharder -d $TEST -f diehard.numbers -g 202 >> results.diehard
  done
  echo ""
fi
