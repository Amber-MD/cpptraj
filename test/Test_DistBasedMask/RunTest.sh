#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in First.pdb.1 Second.pdb.1 Third.pdb.1 Fourth.pdb.1

TESTNAME='Distance-based atom mask tests'
Requires netcdf

INPUT='-i mask.in'

cat > mask.in <<EOF
parm ../tz2.parm7
reference ../tz2.nc 1 [FIRST]
reference ../tz2.nc lastframe [LAST]
trajin ../tz2.nc 1 1

activeref ref [FIRST]
# Atoms outside 5 Ang from residue 2
mask :2>@5.0 maskpdb First.pdb
# Atoms within 5 Ang of residue 2
mask :2<@5.0 maskpdb Second.pdb

activeref ref [LAST]
# Residues within 5 Ang of residue 2
strip !(:2<:5.0)
outtraj Third.pdb.1 pdb noter
unstrip
# Residues outside 5 Ang from residue 2
strip !(:2>:5.0)
outtraj Fourth.pdb.1 pdb noter

run
EOF
RunCpptraj "$TESTNAME"
DoTest First.pdb.1.save First.pdb.1
DoTest Second.pdb.1.save Second.pdb.1
DoTest Third.pdb.1.save Third.pdb.1
DoTest Fourth.pdb.1.save Fourth.pdb.1

EndTest

exit 0
