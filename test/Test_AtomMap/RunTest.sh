#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles atommap.in initial.mol2 atommap.dat reordered.pdb reordered.mol2 fit.mol2 rmsd.dat 

# Test 1
cat > atommap.in <<EOF
noprogress
parm xtallig.mol2
reference xtallig.mol2

parm start.mol2
reference start.mol2 parmindex 1

outtraj initial.mol2 mol2

# xtallig -> start
atommap xtallig.mol2 start.mol2 mapout atommap.dat

trajin xtallig.mol2
outtraj reordered.pdb pdb chainid X model
outtraj reordered.mol2 mol2
rms 1g9v refindex 1 out rmsd.dat
trajout fit.mol2 mol2
EOF
INPUT="-i atommap.in"
RunCpptraj "AtomMap Test"

DoTest initial.mol2.save initial.mol2
DoTest atommap.dat.save atommap.dat
DoTest reordered.pdb.save reordered.pdb
DoTest reordered.mol2.save reordered.mol2
DoTest rmsd.dat.save rmsd.dat
DoTest fit.mol2.save fit.mol2

cat > atommap.in <<EOF
parm xtallig.mol2
reference xtallig.mol2

parm start.mol2
reference start.mol2 parm start.mol2

atommap xtallig.mol2 start.mol2 rmsfit rmsout rmsout.dat 1g9v

trajin xtallig.mol2
EOF
RunCpptraj "AtomMap with 'rmsfit'"
DoTest rmsd.dat.save rmsout.dat

EndTest

exit 0
