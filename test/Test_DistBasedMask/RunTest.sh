#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in First.dat First.pdb.1 Second.dat Second.pdb.1 \
  Third.dat Third.pdb.1 Fourth.pdb.1

TESTNAME='Distance-based atom mask tests'
Requires netcdf

INPUT='-i mask.in'

cat > mask.in <<EOF
parm ../tz2.parm7
reference ../tz2.nc 1 [FIRST]
trajin ../tz2.nc 1 1

# Atoms outside 5 Ang from residue 2
mask :2>@5.0 maskout First.dat maskpdb First.pdb
# Atoms within 5 Ang of residue 2
mask :2<@5.0 maskout Second.dat maskpdb Second.pdb
run
EOF
RunCpptraj "$TESTNAME, by atom."
DoTest First.pdb.1.save First.pdb.1
DoTest Second.pdb.1.save Second.pdb.1

cat > mask.in <<EOF
parm ../tz2.parm7
reference ../tz2.nc lastframe [LAST]
trajin ../tz2.nc 1 1

# Residues within 5 Ang of residue 2
strip !(:2<:5.0)
outtraj Third.pdb.1 pdb noter
unstrip
# Residues outside 5 Ang from residue 2
strip !(:2>:5.0)
outtraj Fourth.pdb.1 pdb noter
run
EOF
RunCpptraj "$TESTNAME, by residue."
DoTest Third.pdb.1.save Third.pdb.1
DoTest Fourth.pdb.1.save Fourth.pdb.1

EndTest

exit 0
