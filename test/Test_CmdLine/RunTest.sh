#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cmd.in test1.crd.save test1.crd

TESTNAME='Command line tests'
# Required environment 
Requires notparallel

INPUT='-i cmd.in'

# First test without command line
cat > cmd.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd 2 13 3
trajin ../tz2.crd 98 last
trajout test1.crd.save
EOF
RunCpptraj "Command line test, part 1"

# Now test using the command line. Hijack the INPUT variable.
INPUT="-p ../tz2.parm7 -y ../tz2.crd -ya \"2 13 3\" -y ../tz2.crd -ya \"98 last\" -x test1.crd"
RunCpptraj "Command line test, part 2"
DoTest test1.crd.save test1.crd


EndTest

exit 0
