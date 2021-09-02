#!/bin/bash

make clean
make

MODE=2
TESTLIST='0 1 2 3 4'

for TEST in $TESTLIST ; do
  echo "Test rng $TEST"
  ./a.out --mode $MODE -r $TEST > results.$MODE.$TEST
done
