#!/bin/bash

if [ ! -f 'UnitMaster.sh' ] ; then
  echo "Fatal Error: UnitMaster.sh not present" > /dev/stderr
  exit 1
fi
CPPTRAJ_TEST_MODE='master'
. UnitMaster.sh
exit 0
