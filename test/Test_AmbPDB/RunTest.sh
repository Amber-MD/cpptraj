#!/bin/bash

. ../MasterTest.sh

CleanFiles out.pdb out.mol2

if [ ! -f "$AMBPDB" ] ; then
  echo "Warning: AMBPDB is not correctly set. Skipping."
  if [[ ! -z $TEST_RESULTS ]] ; then
    echo "Warning: AMBPDB is not correctly set. Skipping." >> $TEST_RESULTS
  fi
  EndTest
  exit 0
else
  $VALGRIND $AMBPDB -p ../tz2.parm7 -c ../tz2.rst7 > out.pdb 2> $ERROR
  DoTest out.pdb.save out.pdb

  $VALGRIND $AMBPDB -p ../tz2.parm7 < ../tz2.rst7 > out.pdb 2>> $ERROR
  DoTest out.pdb.save out.pdb
  if [[ ! -z $AMBERHOME ]] ; then
    $VALGRIND $AMBPDB -p ../tz2.parm7 -mol2 -sybyl < ../tz2.rst7 > out.mol2 2>> $ERROR
    DoTest ../Test_Mol2/test2.mol2.save out.mol2
  else
    echo "Amber to SYBYL atom type conversion test requires AMBERHOME be set."
  fi
fi

EndTest
exit 0
