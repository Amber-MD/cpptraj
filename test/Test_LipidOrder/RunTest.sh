#!/bin/bash

. ../MasterTest.sh

in=order.in
out1=sn1_est.dat
oute1=e2e_sn1.dat
out2=sn2_est.dat
oute2=e2e_sn2.dat
out3=sn1_dir.dat
out4=sn2_dir.dat

CleanFiles $in $out1 $oute1 $out2 $oute2 $out3 $out4

NotParallel "Lipid Order Parameter Test."
if [[ $? -ne 0 ]] ; then
  echo ""
  exit 0
fi

INPUT="-i $in"


del='delta 0.1'

cat > $in <<EOF
# crd/top courtesy of Callum Dickson, Imperial College London
parm ../DOPC.parm7
trajin ../DOPC.rst7

lipidorder out $out1 z unsat ":OL@C19 | :OL@C110" \
 taildist $oute1 $del tailstart ":OL@C12" tailend ":OL@C118" \
 ":OL@C12" ":OL@C13" ":OL@C14" ":OL@C15" ":OL@C16" ":OL@C17" \
 ":OL@C18" ":OL@C19" ":OL@C110" ":OL@C111" ":OL@C112" ":OL@C113" \
 ":OL@C114" ":OL@C115" ":OL@C116" ":OL@C117" ":OL@C118"

lipidorder out $out2 z unsat ":OL2@C19 | :OL2@C110" \
  taildist $oute2 $del tailstart ":OL2@C12" tailend ":OL2@C118" \
  ":OL2@C12" ":OL2@C13" ":OL2@C14" ":OL2@C15" ":OL2@C16" ":OL2@C17" \
  ":OL2@C18" ":OL2@C19" ":OL2@C110" ":OL2@C111" ":OL2@C112" ":OL2@C113" \
  ":OL2@C114" ":OL2@C115" ":OL2@C116" ":OL2@C117" ":OL2@C118"

lipidorder out $out3 z scd \
  ":OL@C12" ":OL@H2R" ":OL@H2S" ":OL@C13" ":OL@H3R" ":OL@H3S" \
  ":OL@C14" ":OL@H4R" ":OL@H4S" ":OL@C15" ":OL@H5R" ":OL@H5S" \
  ":OL@C16" ":OL@H6R" ":OL@H6S" ":OL@C17" ":OL@H7R" ":OL@H7S" \
  ":OL@C18" ":OL@H8R" ":OL@H8S" ":OL@C19" ":OL@H9R" ":OL@H9R" \
  ":OL@C110" ":OL@H10R" ":OL@H10R" ":OL@C111" ":OL@H11R" ":OL@H11S" \
  ":OL@C112" ":OL@H12R" ":OL@H12S" ":OL@C113" ":OL@H13R" ":OL@H13S" \
  ":OL@C114" ":OL@H14R" ":OL@H14S" ":OL@C115" ":OL@H15R" ":OL@H15S" \
  ":OL@C116" ":OL@H16R" ":OL@H16S" ":OL@C117" ":OL@H17R" ":OL@H17S" \
  ":OL@C118" ":OL@H18R" ":OL@H18S"

lipidorder out $out4 z scd \
  ":OL2@C12" ":OL2@H2R" ":OL2@H2S" ":OL2@C13" ":OL2@H3R" ":OL2@H3S" \
  ":OL2@C14" ":OL2@H4R" ":OL2@H4S" ":OL2@C15" ":OL2@H5R" ":OL2@H5S" \
  ":OL2@C16" ":OL2@H6R" ":OL2@H6S" ":OL2@C17" ":OL2@H7R" ":OL2@H7S" \
  ":OL2@C18" ":OL2@H8R" ":OL2@H8S" ":OL2@C19" ":OL2@H9R" ":OL2@H9R" \
  ":OL2@C110" ":OL2@H10R" ":OL2@H10R" ":OL2@C111" ":OL2@H11R" ":OL2@H11S" \
  ":OL2@C112" ":OL2@H12R" ":OL2@H12S" ":OL2@C113" ":OL2@H13R" ":OL2@H13S" \
  ":OL2@C114" ":OL2@H14R" ":OL2@H14S" ":OL2@C115" ":OL2@H15R" ":OL2@H15S" \
  ":OL2@C116" ":OL2@H16R" ":OL2@H16S" ":OL2@C117" ":OL2@H17R" ":OL2@H17S" \
  ":OL2@C118" ":OL2@H18R" ":OL2@H18S"
EOF

RunCpptraj "Lipid Order Parameter Test."
DoTest ${out1}.save $out1
DoTest ${oute1}.save $oute1
DoTest ${out2}.save $out2
DoTest ${oute2}.save $oute2
DoTest ${out3}.save $out3
DoTest ${out4}.save $out4
CheckTest
EndTest

exit 0
