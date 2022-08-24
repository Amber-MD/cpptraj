#!/bin/bash

. ../MasterTest.sh

# NOTE: This test is a little different in that it does not call RunCpptraj
#       directly.
#       The ndiff functionality returns number of differences as an exit
#       status.

CleanFiles diff.out

TESTNAME='Test ndiff functionality'

((PROGCOUNT++))

echo ""
echo "  CPPTRAJ: $TESTNAME"
if [ -z "$CPPTRAJ_DACDIF" ] ; then
  OUT "  CPPTRAJ: $TESTNAME"
fi

$CPPTRAJ --ndiff -v ABSERR=0.000001 ../Test_RMSD/vecs.dat.save vecs.dat > diff.out

# There should be 9 differences.
STATUS=$?
if [ $STATUS -ne 9 ] ; then
  ProgramError "cpptraj exited with status $STATUS"
fi

DoTest diff.out.save diff.out

EndTest
