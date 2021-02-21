#!/bin/bash

COUNT=100000000

make rngtest

# NOTE: Omitting test 200 for now. Gives following message:
#Error:  Can only test distribution of positive ntuples.
#        Use -n ntuple for 0 < ntuple.
#        Read test description with dieharder -d 200 -h.

TESTNUMS='0 1 2 3 4 8 9 10 11 12 13 15 16 17 100 101 102 201 202 203 204 205 206 207 208 209'

TYPELIST=''
DIEHARD='no'
while [ ! -z "$1" ] ; do
  case "$1" in
    '-type' ) shift; TYPELIST="$TYPELIST $1" ;;
  esac
  shift
done
if [ -z "$TYPELIST" ] ; then
  TYPELIST='0 1 2 3 4'
  DIEHARD='yes'
fi
echo "Type list: $TYPELIST"

for CPPTRAJTYPE in $TYPELIST ; do
  OUTFILE=results.$CPPTRAJTYPE
  echo "Cpptraj type $CPPTRAJTYPE"
  echo "Cpptraj type $CPPTRAJTYPE" > $OUTFILE
  TMPFILE=temp.$CPPTRAJTYPE
  ./rngtest -t $COUNT -S 1 -r $CPPTRAJTYPE > $TMPFILE
  for TEST in $TESTNUMS ; do
    echo "Test $TEST"
    dieharder -d $TEST -f $TMPFILE -g 202 >> $OUTFILE
  done
  echo ""
done

if [ "$DIEHARD" = 'yes' ] ; then
  echo "Diehard numbers"
  echo "Diehard numbers" > results.diehard
  dieharder -S 1 -B -o -t $COUNT > diehard.numbers
  for TEST in $TESTNUMS ; do
    echo "Test $TEST"
    dieharder -d $TEST -f diehard.numbers -g 202 >> results.diehard
  done
  echo ""
fi
