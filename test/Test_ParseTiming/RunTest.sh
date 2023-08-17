#!/bin/bash

. ../MasterTest.sh

TESTNAME='Parse timing tests'

CleanFiles parse.in parse.dat parse.out

INPUT='-i parse.in'

UNITNAME='Parse timing test'
cat > parse.in <<EOF
parsetiming hpod*.out sortby time out parse.dat includebad groupout parse.out grouptype kind showdetails
list
EOF
RunCpptraj "$UNITNAME"
DoTest parse.dat.save parse.dat
DoTest parse.out.save parse.out

EndTest

