#!/bin/bash

. ../MasterTest.sh

CleanFiles create.in crd?.dat crd?.gnu crd?.rst7 crd.rst7

INPUT="-i create.in"

# Test COORDS creation and processing
cat > create.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
atomicfluct out crd1.dat byatom bfactor
createcrd crd1
run
crdfluct crdset crd1 out crd1.dat window 20 bfactor
runanalysis
EOF
RunCpptraj "COORDS data set creation and CRDFLUCT test."
DoTest crd1.dat.save crd1.dat
 
# Test velocities
cat > create.in <<EOF
#parm trpzip2.ff10.mbondi.parm7
#trajin tz2.220.rst7
parm ../rGACC.full.parm7
trajin ../rGACC.full.nc
trajout crd1.rst7 
trajout crd0.rst7 title " " nobox
createcrd crd1
run
crdout crd1 crd.rst7 title " " nobox time0 50961
EOF
RunCpptraj "COORDS data set creation with velocities test."
DoTest crd0.rst7 crd.rst7
DoTest crd1.rst7.save crd1.rst7

# Test TRAJ creation and processing
cat > create.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
loadtraj name crd1
atomicfluct out crd1.dat byatom bfactor
run
crdfluct crdset crd1 out crd1.dat window 20 bfactor
runanalysis
EOF
RunCpptraj "TRAJ data set creation and CRDFLUCT test."
DoTest crd1.dat.save crd1.dat

# Test TRAJ creation with velocities
cat > create.in <<EOF
parm ../rGACC.full.parm7
loadtraj ../rGACC.full.nc name crd2
crdout crd2 crd2.rst7 time0 50961
EOF
RunCpptraj "TRAJ data set creation with velocities test."
DoTest crd1.rst7.save crd2.rst7

Disable() {
cat > create.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
atomicfluct out crd1.dat byatom bfactor
createcrd crd1
crdfluct crd1 out crd1.gnu window 20
#list
run
#gnuplot crd1.gnu
#list analysis
clear analysis
crdfluct crd1 out crd2.gnu window 5
runanalysis
list datafile
EOF
#RunCpptraj

cat > create.in <<EOF
parm ../dna30.parm7 
trajin ../dna30.fs11.0.nc 
parminfo  
molinfo !:WAT
strip :WAT,@Na+  
select :WAT,@Na+
list actions
clear actions
strip :WAT
strip @Na+
select @Na+
list actions
createcrd dna
run
list datasets
help nastruct
reference ../dna30.fs11.0.nc 1 
list
help nastruct
crdaction dna nastruct reference
EOF
#RunCpptraj "DNA"

cat > create.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
#strip :WAT outprefix nowat
createcrd crd1
run
crdout crd1 tz2.nowat.rst7 restart onlyframes 1
crdaction crd1 autoimage
crdout crd1 tz2.nowat.nc netcdf
EOF
#RunCpptraj "CRD with box"
}
EndTest

exit 0
