#!/bin/bash

. ../MasterTest.sh

CleanFiles mdXtal.pdb xtals.in test.out

TESTNAME='XtalSymm tests'
Requires netcdf
INPUT="-i xtals.in"

# Simple reimaging of just the asymmetric unit
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1
trajout mdXtal.pdb
EOF
RunCpptraj "XtalSymm Reference"
DoTest mdAsuOnly.pdb.save mdXtal.pdb

# Reimaging of the entire solvent, atom by atom
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1 collect
trajout mdXtal.pdb
EOF
RunCpptraj "XtalSymm Reimaging by Atom"
DoTest mdSolventByAtom.pdb.save mdXtal.pdb

# Reimaging of the entire solvent, atom by atom
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1 collect centroid
trajout mdXtal.pdb
EOF
RunCpptraj "XtalSymm Reimaging by Molecule"
DoTest mdSolventByMolecule.pdb.save mdXtal.pdb

EndTest

exit 0
