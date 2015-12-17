#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles dummy.rst7 dummy.pdb strip.in *.tz2.truncoct.parm7 Complex.crd Receptor.crd Ligand.crd

NotParallel
# Strip Test
cat > strip.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.crd 1 1 parm ../tz2.truncoct.parm7
strip :14-16,18-99999 outprefix strip
trajout dummy.pdb pdb parm ../tz2.truncoct.parm7 chainid X nobox
trajout dummy.rst7 restart parm ../tz2.truncoct.parm7
EOF
INPUT="-i strip.in"
RunCpptraj "One frame strip command test."
DoTest dummy.pdb.save dummy.pdb
DoTest dummy.rst7.save dummy.rst7
# Tell diff to ignore the VERSION line
DoTest strip.tz2.truncoct.parm7.save strip.tz2.truncoct.parm7 -I %VERSION 
CheckTest

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
INPUT="-i strip.in"
RunCpptraj "Unstrip (Lig/Rec/Complex) command test."
# Tell diff to ignore the VERSION line
DoTest strip.tz2.truncoct.parm7.save complex.tz2.truncoct.parm7 -I %VERSION
DoTest receptor.tz2.truncoct.parm7.save receptor.tz2.truncoct.parm7 -I %VERSION
DoTest ligand.tz2.truncoct.parm7.save ligand.tz2.truncoct.parm7 -I %VERSION
DoTest Ligand.crd.save Ligand.crd
DoTest Receptor.crd.save Receptor.crd
DoTest Complex.crd.save Complex.crd
CheckTest

EndTest

exit 0
