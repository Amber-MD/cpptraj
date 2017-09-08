#!/bin/bash

if [ ! -f 'MasterTest.sh' ] ; then
  echo "Fatal Error: MasterTest.sh not present" > /dev/stderr
  exit 1
fi
CPPTRAJ_TEST_MODE='master'
. MasterTest.sh
exit 0
