#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cmd.in test1.crd.save test1.crd res res.save

TESTNAME='Command line tests'
# Required environment 
Requires maxthreads 8


# First test without command line
cat > cmd.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd 2 13 3
trajin ../tz2.crd 98 last
trajout test1.crd.save title MyTitle
EOF
INPUT='-i cmd.in'
RunCpptraj "Command line test, part 1"

# Now test using the command line. Hijack the INPUT variable.
INPUT="-p ../tz2.parm7 -y ../tz2.crd -ya \"2 13 3\" -y ../tz2.crd -ya \"98 last\" -x test1.crd -xa \"title MyTitle\""
RunCpptraj "Command line test, part 2"
DoTest test1.crd.save test1.crd

# Test loading multiple topology files
cat > cmd.in <<EOF
parm ../tz2.parm7
parm ../tz2.pdb
resinfo out res.save parmindex 0
resinfo out res.save parmindex 1
EOF
INPUT='-i cmd.in'
RunCpptraj "Command line wildcard test, part 1"
cat > cmd.in <<EOF
resinfo out res parmindex 0
resinfo out res parmindex 1
EOF
INPUT='../tz2.p* cmd.in'
RunCpptraj "Command line wildcard test, part 2"
DoTest res.save res

EndTest

exit 0
