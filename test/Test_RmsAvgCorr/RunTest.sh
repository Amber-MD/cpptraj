#!/bin/bash

. ../MasterTest.sh

CleanFiles rms.in rmscorr.dat rmscorr.10.dat rmscorr.first.dat
CheckNetcdf
INPUT="-i rms.in"

# Single reference
cat > rms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
strip !(:2-12@CA)
parm strip.tz2.parm7 [strip]
rmsavgcorr out rmscorr.dat R2-12 reference avg.CA.rst7 parm [strip]
rmsavgcorr out rmscorr.10.dat offset 10 reference avg.CA.rst7 parm [strip]
EOF
RunCpptraj "RmsAvgCorr Reference"
DoTest rmscorr.dat.save rmscorr.dat
DoTest rmscorr.10.dat.save rmscorr.10.dat

# First running avgd frame
cat > rms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
strip !(:2-12@CA)
rmsavgcorr out rmscorr.first.dat First first
EOF
RunCpptraj "RmsAvgCorr First"
DoTest rmscorr.first.dat.save rmscorr.first.dat
CheckTest
EndTest

exit 0
