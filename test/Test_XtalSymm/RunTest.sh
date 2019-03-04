#!/bin/bash

. ../MasterTest.sh

CleanFiles xtals.in mdAsuOnly.pdb mdSolventByAtom.pdb mdSolventByMolecule.pdb

TESTNAME='XtalSymm tests'
Requires netcdf notparallel
INPUT="-i xtals.in"

# Simple reimaging of just the asymmetric unit
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1
trajout mdAsuOnly.pdb
EOF
RunCpptraj "XtalSymm Reference"
DoTest mdAsuOnly.pdb.save mdAsuOnly.pdb

# Reimaging of the entire solvent, atom by atom
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1 collect
trajout mdSolventByAtom.pdb
EOF
RunCpptraj "XtalSymm Reimaging by Atom"
DoTest mdSolventByAtom.pdb.save mdSolventByAtom.pdb

# Reimaging of the entire solvent, atom by atom
cat > xtals.in << EOF
parm ../x6dky.parm7
trajin ../mdXtal.nc
reference ../mdXtal.inpcrd
xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1 collect centroid
trajout mdSolventByMolecule.pdb
EOF
RunCpptraj "XtalSymm Reimaging by Molecule"
DoTest mdSolventByMolecule.pdb.save mdSolventByMolecule.pdb

EndTest

exit 0
