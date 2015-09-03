#!/bin/bash

. ../MasterTest.sh

CleanFiles out.pdb

if [ ! -f "$AMBPDB" ] ; then
  echo "AMBPDB is not correctly set."
  exit 1
else
  $VALGRIND $AMBPDB -p ../tz2.parm7 -c ../tz2.rst7 > out.pdb 2> $ERROR
  DoTest out.pdb.save out.pdb

  $VALGRIND $AMBPDB -p ../tz2.parm7 < ../tz2.rst7 > out.pdb 2>> $ERROR
  DoTest out.pdb.save out.pdb
fi

EndTest
exit 0
