#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles mask.in First.pdb.1 Second.pdb.1 Third.pdb.1 Fourth.pdb.1 \
           Fifth.pdb.1 Sixth.pdb.1

# NOTE: This also tests activeref
TESTNAME='Distance-based atom mask tests'
Requires netcdf maxthreads 1

INPUT='-i mask.in'

cat > mask.in <<EOF
parm ../tz2.parm7
reference ../tz2.nc 1 [FIRST]
reference ../tz2.nc lastframe [LAST]
trajin ../tz2.nc 1 1

activeref ref [FIRST]
# Atoms outside 5 Ang from residue 2
mask :2>@5.0 maskpdb First.pdb trajargs "chainid ' '"
# Atoms within 5 Ang of residue 2
mask :2<@5.0 maskpdb Second.pdb trajargs "chainid ' '"

activeref ref [LAST]
# Residues within 5 Ang of residue 2
strip !(:2<:5.0)
outtraj Third.pdb.1 pdb noter chainid ' '
unstrip
# Residues outside 5 Ang from residue 2
strip !(:2>:5.0)
outtraj Fourth.pdb.1 pdb noter chainid ' '

run
EOF
RunCpptraj "$TESTNAME (atom, residue)"
DoTest First.pdb.1.save First.pdb.1
DoTest Second.pdb.1.save Second.pdb.1
DoTest Third.pdb.1.save Third.pdb.1
DoTest Fourth.pdb.1.save Fourth.pdb.1

cat > mask.in <<EOF
parm ../DOPC.parm7
trajin ../DOPC.rst7
reference ../DOPC.rst7

# Molecules within 3 Ang of residue 1
strip !(:1<^3.0)
outtraj Fifth.pdb.1 pdb noter chainid ' '
unstrip
# Molecules outside 10 Ang of residue 24
strip (!(:24>^22.0))|:WAT
outtraj Sixth.pdb.1 pdb noter chainid ' '

run
EOF
RunCpptraj "$TESTNAME (molecule)"
DoTest Fifth.pdb.1.save Fifth.pdb.1
DoTest Sixth.pdb.1.save Sixth.pdb.1

EndTest

exit 0
