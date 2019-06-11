#!/bin/bash

. ../MasterTest.sh

CleanFiles for.in TRP.vec.dat TRP.rms.dat TRP.CA.dist.dat TRP.tocenter.dat \
           nh.dat rms.nofit.dat last10.dat distance.dat

TESTNAME='Loop tests'
Requires netcdf maxthreads 10

INPUT='-i for.in'
cat > for.in <<EOF
set FNAME= ../tz2.nc
parm ../tz2.parm7
trajin \$FNAME 1 10
for residues T inmask :TRP
  vector \$T center out TRP.vec.dat
  rms first \$T out TRP.rms.dat
  for atoms A0 inmask @CA
    distance \$T \$A0 out TRP.CA.dist.dat
  done
  distance \$T :1-12 out TRP.tocenter.dat
done

for atoms A0 inmask :2-4@N atoms A1 inmask :2-4@H
  distance \$A0 \$A1 out nh.dat
done

rms :1-12&!@H=
for residues R1 inmask :1-12 r=1;r++
  rms R\$r \$R1&!@H= nofit out rms.nofit.dat
done
show

# Print info for the last 10 atoms. This tests using data set values
# as script variables and replacement of multiple script variables
# in a single argument.
set Natom = atoms inmask *
last10 = \$Natom - 10
show
atoms "@\$last10 - \$Natom" out last10.dat
run

EOF
RunCpptraj "$TESTNAME"
DoTest TRP.vec.dat.save TRP.vec.dat
DoTest TRP.rms.dat.save TRP.rms.dat
DoTest TRP.CA.dist.dat.save TRP.CA.dist.dat
DoTest TRP.tocenter.dat.save TRP.tocenter.dat
DoTest nh.dat.save nh.dat
DoTest rms.nofit.dat.save rms.nofit.dat
DoTest last10.dat.save last10.dat

UNITNAME='Test reading comma-separated list of strings'
CheckFor maxthreads 4
if [ $? -eq 0 ] ; then
  cat > for.in <<EOF
parm ../tz2.parm7
# Test reading in a comma-separated list of strings
clear trajin
for TRAJ in ../tz2.nc,../tz2.crd,../tz2.pdb,../tz2.rst7
  trajin \$TRAJ lastframe
done
distance D1-12 :1 :12 out distance.dat
run
EOF
  RunCpptraj "$UNITNAME"
  DoTest distance.dat.save distance.dat
fi

EndTest
exit 0
