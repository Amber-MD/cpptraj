#!/bin/bash

. ../MasterTest.sh

CleanFiles nc.in nc.hp1.ca.dat nc.hp2.ca.dat nc.all.res.dat cmap.dat \
  native.cmap.gnu nonnative.cmap.gnu native.resmap.gnu nonnative.resmap.gnu \
  nc1.pdb nc2.contacts.dat nc2.res.dat NC2.series.dat

INPUT="-i nc.in"
cat > nc.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc

reference ../DPDP.nc 1

nativecontacts name NC1 :1-12@CA out nc.hp1.ca.dat mindist maxdist map mapout cmap.gnu contactpdb nc1.pdb
nativecontacts name NC2 :10-14@CA :17-21@CA out nc.hp2.ca.dat mindist maxdist writecontacts nc2.contacts.dat resout nc2.res.dat series seriesout NC2.series.dat
nativecontacts name NC3 :1-21&!@H= byresidue out nc.all.res.dat mindist maxdist distance 3.0 reference map mapout resmap.gnu
EOF
RunCpptraj "NativeContacts test."
DoTest nc.hp1.ca.dat.save nc.hp1.ca.dat
DoTest nc1.pdb.save nc1.pdb
DoTest nc.hp2.ca.dat.save nc.hp2.ca.dat
DoTest nc2.contacts.dat.save nc2.contacts.dat
DoTest nc2.res.dat.save nc2.res.dat
DoTest NC2.series.dat.save NC2.series.dat
DoTest nc.all.res.dat.save nc.all.res.dat
DoTest native.resmap.gnu.save native.resmap.gnu
DoTest nonnative.resmap.gnu.save nonnative.resmap.gnu
EndTest
exit 0 
