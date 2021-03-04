#!/bin/bash

# Assuming DO_PARALLEL is set with the command to execute mpirun,
# determine how many processes.

if [ -z "`which grep`" ] ; then
  >&2 echo "Error: nproc.sh relies on grep"
  exit 1
fi
if [ -z "`which wc`" ] ; then
  >&2 echo "Error: nproc.sh relies on wc"
  exit 1
fi

NPROC=0

if [ ! -z "$DO_PARALLEL" ] ; then
  NPROC=`$DO_PARALLEL echo "cpptraj_nproc_test" | grep "cpptraj_nproc_test" | wc -l`
fi
echo $NPROC
exit 0
