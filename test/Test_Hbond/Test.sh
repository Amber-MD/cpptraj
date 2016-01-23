#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat solvhb.dat solvavg.dat \
           nbb.dat hbavg.dat solutehb.agr lifehb.gnu avg.lifehb.gnu max.lifehb.gnu \
           crv.lifehb.gnu hb?.dat hbond.mol.dat mol.avg.dat solutehb.dat

INPUT="-i hbond.in"
CheckNetcdf

# Solute-solute, series output, lifetime analysis
cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc

hbond BB1 donormask :14@N acceptormask :17@O series 

hbond BB2 donormask :13@N acceptormask :3@O series

run
write solutehb.dat BB*[solutehb]
EOF
RunCpptraj "Solute Hbond test."
DoTest solutehb.dat.save solutehb.dat

EndTest

exit 0
