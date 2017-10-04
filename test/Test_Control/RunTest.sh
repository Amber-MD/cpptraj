#!/bin/bash

. ../MasterTest.sh

CleanFiles for.in TRP.vec.dat TRP.rms.dat TRP.CA.dist.dat TRP.tocenter.dat \
           nh.dat rms.nofit.dat

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
EOF
RunCpptraj "Loop tests"
DoTest TRP.vec.dat.save TRP.vec.dat
DoTest TRP.rms.dat.save TRP.rms.dat
DoTest TRP.CA.dist.dat.save TRP.CA.dist.dat
DoTest TRP.tocenter.dat.save TRP.tocenter.dat
DoTest nh.dat.save nh.dat
DoTest rms.nofit.dat.save rms.nofit.dat

#DoTest test.out.save test.out

EndTest
exit 0
