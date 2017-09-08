#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles dummy.rst7 dummy.pdb.1 strip.in *.tz2.truncoct.parm7 Complex.crd \
           Receptor.crd Ligand.crd res1.tz2.crd
INPUT="-i strip.in"

# NOTE: strip is also tested in Test_Center
UNITNAME='One frame strip command test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > strip.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.crd 1 1 parm ../tz2.truncoct.parm7
strip :14-16,18-99999 outprefix strip
trajout dummy.pdb pdb multi parm ../tz2.truncoct.parm7 chainid X nobox
trajout dummy.rst7 restart parm ../tz2.truncoct.parm7
EOF
  RunCpptraj "$UNITNAME"
  DoTest dummy.pdb.save dummy.pdb.1
  DoTest dummy.rst7.save dummy.rst7
  # Tell diff to ignore the VERSION line
  DoTest strip.tz2.truncoct.parm7.save strip.tz2.truncoct.parm7 -I %VERSION
fi

UNITNAME='Unstrip (Lig/Rec/Complex) command test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  # Unstrip Test
  cat > strip.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.crd 1 1 parm ../tz2.truncoct.parm7

# Count res 17 as ligand, 1-13 as receptor, 1-13,17 as complex

# Complex
strip !(:1-13,17) outprefix complex
outtraj Complex.crd
unstrip

# Receptor
strip :WAT outprefix receptor
outtraj Receptor.crd
unstrip

# Ligand
strip !(:17) outprefix ligand
outtraj Ligand.crd
EOF
  RunCpptraj "$UNITNAME"
  # Tell diff to ignore the VERSION line
  DoTest strip.tz2.truncoct.parm7.save complex.tz2.truncoct.parm7 -I %VERSION
  DoTest receptor.tz2.truncoct.parm7.save receptor.tz2.truncoct.parm7 -I %VERSION
  DoTest ligand.tz2.truncoct.parm7.save ligand.tz2.truncoct.parm7 -I %VERSION
  DoTest Ligand.crd.save Ligand.crd
  DoTest Receptor.crd.save Receptor.crd
  DoTest Complex.crd.save Complex.crd
fi

# Strip test that will work in parallel
UNITNAME='Multiple frame strip command test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > strip.in <<EOF
noprogress
parm ../tz2.parm7
trajin ../tz2.nc
strip !(:1) nobox
trajout res1.tz2.crd
EOF
  RunCpptraj "Multi frame strip command test."
  DoTest res1.tz2.crd.save res1.tz2.crd
fi
EndTest

exit 0
