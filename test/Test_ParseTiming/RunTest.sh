#!/bin/bash

. ../MasterTest.sh

TESTNAME='Parse timing tests'

CleanFiles parse.in parse.dat parse.out

INPUT='-i parse.in'

UNITNAME='Parse timing test'
cat > parse.in <<EOF
parsetiming hpod13-n.16-N.4-c.4.out hpod13-n.4-N.1-c.4.out hpod13-n.8-N.2-c.4.out \
  sortby time out parse.dat includebad groupout parse.out grouptype kind showdetails name Run
list
EOF
RunCpptraj "$UNITNAME"
DoTest parse.dat.save parse.dat
DoTest parse.out.save parse.out

EndTest

